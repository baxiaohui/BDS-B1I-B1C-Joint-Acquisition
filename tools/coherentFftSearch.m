function [maps,best,meta] = coherentFftSearch(parts,fs,fdC,cfg)
% Exact partitioned matched filtering. Complex correlations are added BEFORE
% powers/normalization, including the final fractional-period partition.
% refCF has ten integer-ms slot templates. Every slot searches fractional
% receive lag 0:1 ms. RJ is an exact slice of the full C result, not a second
% integration or an independently normalized statistic.
cfg=resolveJointAcquisitionOptions(cfg);
NI=round(fs/1000); NC=10*NI; F=numel(fdC); fdC=fdC(:);
gpu=cfg.gpuEnabled;
if gpu, batch=cfg.frequencyBatchSize;
else, batch=min(2,cfg.frequencyBatchSize); end
bi=ceil(NI/cfg.maxSurfaceCodeColumns); bc=ceil(NC/cfg.maxSurfaceCodeColumns);
bf=ceil(F/cfg.maxSurfaceFrequencyRows);
maps.surfaceI=zeros(ceil(F/bf),ceil(NI/bi),'single');
maps.surfaceC=zeros(ceil(F/bf),ceil(NC/bc),'single');
maps.surfaceJ=maps.surfaceI;
maps.phaseAxisI=centers(NI,bi)/fs*1000;
maps.phaseAxisC=centers(NC,bc)/fs*1000;
maps.fdAxisC=meanBins(fdC,bf); maps.fdAxisI=maps.fdAxisC*1561.098/1575.42;
best.I=struct('score',-Inf,'sample',nan,'frequency',nan);
best.C=best.I; best.J=best.I;
% A frequency-INDEPENDENT background avoids dividing away the coherent
% frequency response of strong signals. Two off-band Doppler probes supply
% empirical correlation-noise levels, without including them as detections.
probeOffset=min(fs/16,max(1000,2/(min(cfg.b1iCoherentMs,cfg.b1cCoherentMs)/1000)));
[noiseI,~,noiseC]=correlate(parts,NI,NC,[fdC(1)-probeOffset fdC(end)+probeOffset],cfg);
noiseI=max(double(host(median(abs(noiseI(:)).^2))),double(eps('single')));
noiseC=max(double(host(median(abs(noiseC(:)).^2))),double(eps('single')));
noiseJ=noiseC; % Common background for the SAME C correlations.
timer=tic; nextUpdate=0;
for first=1:batch:F
    indices=first:min(F,first+batch-1); f=fdC(indices).'; nf=numel(indices);
    [RI,RJ,RC]=correlate(parts,NI,NC,f,cfg);
    power=abs(RI).^2; ratioI=power/noiseI;
    power=abs(RJ).^2; ratioJ=power/noiseJ;
    power=abs(RC).^2; ratioC=power/noiseC;
    if strcmp(cfg.jointPhaseMode,'grid') && all(cfg.jointWeights>0)
        [~,residual]=jointPhaseHypothesis(angle(RJ)-angle(RI), ...
            cfg.jointPhaseMode,cfg.jointPhaseBins);
        residual=cast(residual,'like',ratioI);
        joint=weightedJointStatistic(ratioI,ratioJ,cfg.jointWeights,cfg.jointWeightNormalization,residual);
    else
        joint=weightedJointStatistic(ratioI,ratioJ,cfg.jointWeights,cfg.jointWeightNormalization);
    end
    [value,at]=max(ratioI(:)); [q,j]=ind2sub(size(ratioI),double(host(at)));
    if host(value)>best.I.score
        best.I=struct('score',double(host(value)),'sample',q-1,'frequency',indices(j), ...
            'complexPrompt',double(host(RI(q,j))));
    end
    [value,at]=max(ratioC(:)); [q,j]=ind2sub(size(ratioC),double(host(at)));
    if host(value)>best.C.score
        best.C=struct('score',double(host(value)),'sample',q-1,'frequency',indices(j), ...
            'windowSample',mod(q-1,NI), ...
            'complexPrompt',double(host(RC(q,j))));
    end
    [value,at]=max(joint(:)); [q,j]=ind2sub(size(joint),double(host(at)));
    if host(value)>best.J.score
        theta=nan; difference=nan; phaseError=nan; phaseIndex=nan;
        if all(cfg.jointWeights>0)
            difference=double(host(angle(RJ(q,j))-angle(RI(q,j))));
            [theta,phaseError,phaseIndex]=jointPhaseHypothesis(difference, ...
                cfg.jointPhaseMode,cfg.jointPhaseBins);
        end
        best.J=struct('score',double(host(value)),'sample',q-1,'frequency',indices(j), ...
            'ratioI',double(host(ratioI(q,j))),'ratioC',double(host(ratioJ(q,j))), ...
            'theta',theta,'phaseDifference',difference, ...
            'phaseError',phaseError,'phaseIndex',phaseIndex);
    end
    pooledI=poolPhase(ratioI,bi); pooledC=poolPhase(ratioC,bc); pooledJ=poolPhase(joint,bi);
    for j=1:nf
        row=ceil(indices(j)/bf);
        maps.surfaceI(row,:)=max(maps.surfaceI(row,:),pooledI(:,j).');
        maps.surfaceC(row,:)=max(maps.surfaceC(row,:),pooledC(:,j).');
        maps.surfaceJ(row,:)=max(maps.surfaceJ(row,:),pooledJ(:,j).');
    end
    if toc(timer)>=nextUpdate || indices(end)==F
        fprintf('  coherent search: %d/%d Doppler bins, %.1f s elapsed\n',indices(end),F,toc(timer));
        nextUpdate=toc(timer)+15;
    end
end
meta=struct('fullPhaseSamplesI',NI,'fullPhaseSamplesC',NC,'fullFrequencyBins',F, ...
    'phasePoolingI',bi,'phasePoolingC',bc,'frequencyPooling',bf, ...
    'storage','max pooling for display/storage ONLY; all fine-grid samples searched', ...
    'wallSeconds',toc(timer),'noiseI',noiseI,'noiseCJ',noiseJ,'noiseCF',noiseC, ...
    'backgroundMethod','Fixed empirical off-Doppler probes; no per-frequency renormalization', ...
    'jointWeights',cfg.jointWeights,'jointEffectiveWeights',cfg.jointEffectiveWeights, ...
    'jointWeightNormalization',cfg.jointWeightNormalization, ...
    'jointPhaseMode',cfg.jointPhaseMode,'jointPhaseBins',cfg.jointPhaseEffectiveBins, ...
    'jointPhaseStepDeg',cfg.jointPhaseStepDeg, ...
    'b1cWindowMode','shared-candidate','b1cBackgroundShared',true, ...
    'b1iCoherentMs',cfg.b1iCoherentMs,'b1cCoherentMs',cfg.b1cCoherentMs);
end
function [RI,RJ,RC]=correlate(parts,NI,NC,f,cfg)
gpu=cfg.gpuEnabled; nf=numel(f);
    RI=complex(zeros(NI,nf,'single')); RJ=RI;
    RC=complex(zeros(NC,nf,'single'));
    if gpu, RI=gpuArray(RI); RJ=gpuArray(RJ); RC=gpuArray(RC); end
    for k=1:numel(parts)
        p=parts{k};
        if ischar(p), loaded=load(p,'p'); p=loaded.p; end
        assert(isfield(p,'b1cWindowMode') && strcmp(p.b1cWindowMode,'shared-candidate'), ...
            'Old acquisition partitions: regenerate with prepareCoherentParts.');
        jointSlot=p.cJointSlot;
        if gpu && ~isa(p.xI,'gpuArray')
            for field={'xI','xC','time','refI','refCF'}
                p.(field{1})=gpuArray(p.(field{1}));
            end
        end
        if gpu, fwork=gpuArray(f); else, fwork=f; end
        fi=fwork*(1561.098/1575.42);
        rateI=cfg.dopplerRateHzPerSecond*(1561.098/1575.42);
        % Reduce the phase in double precision before casting to single;
        % the time origin stays at acquisition start for EVERY partition.
        if p.activeI
            phase=single(2*pi*mod(p.time*fi+0.5*rateI*p.time.^2,1));
            y=fft(p.xI.*exp(-1j*phase),p.nfft,1);
            r=ifft(y.*p.refI,[],1); RI=RI+r(1:NI,:);
            clear y r phase;
        end
        if p.activeCJ || p.activeCF
            phase=single(2*pi*mod(p.time*fwork+0.5*cfg.dopplerRateHzPerSecond*p.time.^2,1));
            y=fft(p.xC.*exp(-1j*phase),p.nfft,1);
            if p.activeCF
                for slot=0:9
                    r=ifft(y.*p.refCF(:,slot+1),[],1);
                    % phase = fractional receive lag - primary position.
                    index=mod((0:NI-1)'-slot*NI,NC)+1;
                    RC(index,:)=RC(index,:)+r(1:NI,:);
                end
            end
        end
        clear r y p;
    end
    index=mod((0:NI-1)'-jointSlot*NI,NC)+1;
    RJ=RC(index,:);
end
function out=host(x)
if isa(x,'gpuArray'), out=gather(x); else, out=x; end
end
function pooled=poolPhase(z,b)
n=size(z,1); count=ceil(n/b); missing=count*b-n;
if missing>0, z=[z; -inf(missing,size(z,2),'like',z)]; end
pooled=reshape(max(reshape(z,b,count,[]),[],1),count,[]);
pooled=host(pooled);
end
function x=centers(n,b)
first=1:b:n; last=min(n,first+b-1); x=(first+last)/2-1;
end
function out=meanBins(x,b)
out=zeros(ceil(numel(x)/b),1);
for k=1:numel(out), out(k)=mean(x((k-1)*b+1:min(numel(x),k*b))); end
end
