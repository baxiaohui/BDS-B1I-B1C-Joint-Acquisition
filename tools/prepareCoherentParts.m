function [parts,warp,meta,cleanup] = prepareCoherentParts(cfg,aid,fs)
% Prepare bounded FFT partitions. No capture-length-sized raw LS3W array.
% Known A->B code epochs define a continuous CODE-TIME warp, compensating
% code drift without tracking the weak portion of B.
cfg=resolveJointAcquisitionOptions(cfg);
NI=round(fs/1000); NC=10*NI;
nh=aid.nh.rows; tx=unwrapSeconds(nh(:,1));
if cfg.codeDriftCompensation
    zero=interp1(nh(:,4),tx,cfg.startSecond,'linear');
    warp.u=tx-zero; warp.t=nh(:,4)-cfg.startSecond;
else
    warp.u=[-1;cfg.durationSeconds+1]; warp.t=warp.u;
end
assert(all(isfinite(warp.u)),'Code-time warp lacks the acquisition start.');
first=find(nh(:,4)>=cfg.startSecond,1); boundary=nh(first,4)-cfg.startSecond;
% Every C phase hypothesis uses the SAME fractional-ms receive-window
% start as the corresponding I/joint hypothesis. Search the ten integer
% C slots in the LOCAL replica, rather than moving the samples by 0:9 ms.
countI=sampleCount(warp,boundary,cfg.b1iCoherentMs,fs);
countCJ=sampleCount(warp,boundary,cfg.b1cCoherentMs,fs);
countCF=countCJ;
count=max([countI countCJ countCF]);
maxBlock=max(1,round(cfg.fftSegmentMs*fs/1000));
starts=0:maxBlock:count-1; lens=min(maxBlock,count-starts);
nfft=2.^nextpow2(lens+NI);
estimated=sum(11*8*nfft+24*(lens+NI)); % I + ten C slot templates
limit=cfg.maxWorkingCacheMB*1024^2;
cacheOnGPU=false;
if cfg.gpuEnabled
    device=gpuDevice();
    cacheOnGPU=estimated<min(limit,0.25*device.AvailableMemory);
end
inMemory=estimated<=limit;
work=tempname(cfg.workCacheDir); mkdir(work);
cleanup=onCleanup(@()cleanWork(work,cfg.workCacheDir));
parts=cell(1,numel(starts));
codeI=double(generateB12Icode(cfg.currentPRN)); codeI=codeI(:);
s.codeLength=10230;
dataCode=generateDataBOC11(s,cfg.currentPRN);
pilotCode=generatePilotBOC11(s,cfg.currentPRN);
dataCode=double(dataCode(1:6:end)); dataCode=dataCode(:);
pilotCode=double(pilotCode(1:6:end)); pilotCode=pilotCode(:);
repI=single(codeI(mod(floor((0:NI-1)'*2.046e6/fs),2046)+1));
index=mod(floor((0:NC-1)'*2.046e6/fs),20460)+1;
repD=single(dataCode(index)); repP=single(pilotCode(index));
firstNH=nh(first,:); h=firstNH(5)-1;
txI=firstNH(1);
for k=1:numel(parts)
    offset=starts(k); len=lens(k); sample=(offset+(0:len+NI-1)')./fs;
    time=interp1(warp.u,warp.t,sample,'linear');
    assert(all(isfinite(time)),'Missing code epoch times inside integration window.');
    [xI,xC]=sampleWindow(cfg,time,fs);
    lenI=max(0,min(len,countI-offset));
    lenCJ=max(0,min(len,countCJ-offset));
    lenCF=max(0,min(len,countCF-offset));
    rI=single([]);
    if lenI>0
        n=offset+(0:lenI-1)'; ms=floor(n/NI); allMs=unique(ms);
        signs=lookupTx(aid.dataI,txI+allMs/1000).*lookupTx(aid.nh,txI+allMs/1000);
        rI=repI(mod(n,NI)+1).*single(signs(ms-allMs(1)+1));
    end
    p=struct('xI',xI,'xC',xC,'time',time,'nfft',nfft(k), ...
        'activeI',lenI>0,'activeCJ',lenCJ>0,'activeCF',lenCF>0, ...
        'cJointSlot',mod(h,10),'b1cWindowMode','shared-candidate');
    if cacheOnGPU
        p.xI=gpuArray(xI); p.xC=gpuArray(xC); p.time=gpuArray(time);
        rI=gpuArray(rI);
    end
    p.refI=single([]); p.refCF=single([]);
    if p.activeI, p.refI=conj(fft(rI,nfft(k))); end
    if p.activeCF
        p.refCF=complex(zeros(nfft(k),10,'single'));
        if cacheOnGPU, p.refCF=gpuArray(p.refCF); end
        n=offset+(0:lenCF-1)'; ms=floor(n/NI); allMs=unique(ms);
        signIndex=ms-allMs(1)+1;
        for slot=0:9
            % Shift primary AND its overlay together. The true NH-selected
            % slot has t=txI+(allMs+.5)/1000, exactly the old joint replica.
            t=txI+(allMs+0.5+slot-mod(h,10))/1000;
            nav=lookupTx(aid.dataC,t); secondary=lookupTx(aid.secondary,t);
            index=mod(n+slot*NI,NC)+1;
            rC=b1cMatchedReplica(repD(index),repP(index),nav(signIndex),secondary(signIndex),cfg);
            if cacheOnGPU, rC=gpuArray(rC); end
            p.refCF(:,slot+1)=conj(fft(rC,nfft(k)));
        end
    end
    if inMemory, parts{k}=p;
    else
        parts{k}=fullfile(work,sprintf('part_%06d.mat',k));
        save(parts{k},'p','-v6'); % uncompressed bounded partition, outside output
    end
end
meta=struct('codeSamples',count,'codeSamplesI',countI,'codeSamplesCJ',countCJ, ...
    'codeSamplesCF',countCF,'requestedMsI',cfg.b1iCoherentMs,'requestedMsC',cfg.b1cCoherentMs,'NHIndex',h+1, ...
    'parts',numel(parts),'cacheBytesEstimated',estimated,'cacheOnGPU',cacheOnGPU, ...
    'cacheInMemory',inMemory,'codeDriftCompensation',cfg.codeDriftCompensation, ...
    'b1cBranch',cfg.b1cBranch,'b1cEffectiveDataPilotWeights',cfg.b1cEffectiveDataPilotWeights, ...
    'b1cPilotPhaseDeg',cfg.b1cPilotPhaseDeg,'b1cPilotModel','BOC11 only; no BOC61 matched branch', ...
    'b1cWindowMode','shared-candidate','cJointSlot',mod(h,10));
end
function count=sampleCount(warp,boundary,milliseconds,fs)
u0=interp1(warp.t,warp.u,boundary,'linear');
u1=interp1(warp.t,warp.u,boundary+milliseconds/1000,'linear');
assert(isfinite(u0) && isfinite(u1),'Code-time warp lacks the integration end.');
count=max(1,round((u1-u0)*fs));
end
function out=unwrapSeconds(t)
out=t+604800*cumsum([0;diff(t)<-302400]);
end
function value=lookupTx(aid,t)
x=unwrapSeconds(aid.rows(:,1)); t=t(:);
t=t+604800*round((x(1)-t(1))/604800);
assert(all(t>=x(1)-1e-8 & t<x(end)+aid.step+1e-8),'Missing on-air signs at requested code epochs.');
ix=interp1(x,(1:numel(x))',t+1e-9,'previous','extrap'); value=aid.rows(ix,3);
end
function [i,c]=sampleWindow(cfg,t,fs)
absolute=cfg.fileStartSecond+cfg.startSecond+t;
% Lanczos interpolation, 16 taps, avoids the bandwidth loss of a linear
% fractional-sample shift. Small guard regions cover interpolation support.
begin=max(0,floor(min(absolute)*fs)-9);
last=ceil(max(absolute)*fs)+9;
if isfield(cfg,'weakNoiseModel')
    [ri,rc]=readLabsatWindow(cfg.ls3wFile,begin/fs,(last-begin+1)/fs*1000,cfg.weakNoiseModel);
else
    [ri,rc]=readLabsatWindow(cfg.ls3wFile,begin/fs,(last-begin+1)/fs*1000);
end
position=(absolute-begin/fs)*fs+1;
i=interpolate(ri,position); c=interpolate(rc,position);
end
function y=interpolate(x,position)
base=floor(position); frac=position-base;
y=complex(zeros(size(position),'single')); total=zeros(size(position));
for tap=-7:8
    d=frac-tap; w=sinc(d).*sinc(d/8); w(abs(d)>=8)=0;
    index=base+tap; valid=index>=1 & index<=numel(x);
    y(valid)=y(valid)+single(w(valid)).*single(x(index(valid)));
    total(valid)=total(valid)+w(valid);
end
y=y./single(total);
end
function cleanWork(work,parent)
target=char(java.io.File(work).getCanonicalPath());
prefix=[char(java.io.File(parent).getCanonicalPath()) filesep];
assert(startsWith(lower(target),lower(prefix)),'Unsafe temporary-cache path.');
if ~isfolder(work), return; end
files=dir(fullfile(work,'part_*.mat'));
for k=1:numel(files), delete(fullfile(work,files(k).name)); end
rmdir(work); % nonrecursive; never removes unrelated files
end
