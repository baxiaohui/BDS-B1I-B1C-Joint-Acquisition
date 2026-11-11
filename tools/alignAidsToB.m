function alignment = alignAidsToB(A,B,cfg)
% Match actual AIR bits. No complete B1C/B1I frame is needed in B.
allOutput={};
for signal={'B1I','B1C'}
    s=signal{1}; ar=A.(s).report; br=B.(s).report;
    available=arrayfun(@(r)~isempty(r.files.dataBits),ar);
    prns=intersect([ar(available).PRN],[br.PRN]);
    if ~isempty(cfg.prns), prns=intersect(prns,cfg.prns); end
    assert(~isempty(prns),'No time-synchronized A/B %s satellites.',s);
    records=struct([]); stored=load(B.(s).trackingFile,'trackResults','settings');
    for prn=prns
        a=ar(find([ar.PRN]==prn,1));
        ad=readStrippingAid(a.files.dataBits); ao=readStrippingAid(a.files.overlay);
        t=stored.trackResults(find([stored.trackResults.PRN]==prn,1));
        if strcmp(s,'B1I')
            o=demodB1IAir(t,stored.settings);
            [bd,bo]=matchB1IBits(o,ad,ao,cfg);
            method='AIR_BITS_NH_MATCH';
        else
            model=alignment.B1I([alignment.B1I.PRN]==prn);
            assert(~isempty(model),'B1C C%02d needs same-satellite B1I prefix alignment.',prn);
            [bd,bo]=matchB1C(t,stored.settings,ad,ao,model,cfg);
            method='AIR_BITS_SECONDARY_MATCH';
        end
        [~,ia,ib]=intersect(epochKeys(ao.rows(:,1)),epochKeys(bo.rows(:,1)));
        good=bo.rows(ib,4)>=3 & bo.rows(ib,4)<cfg.B.strongSeconds-0.1;
        ta=ao.rows(ia(good),4); tb=bo.rows(ib(good),4);
        assert(numel(ta)>=100,'C%02d: insufficient common code epochs.',prn);
        [offset,scale,residual,bound]=clockFit(ta,tb,cfg);
        assert(max(abs(residual))*1e6<=cfg.maxAlignmentResidualUs, ...
            '%s C%02d: A/B residual %.2f us exceeds limit.',s,prn,max(abs(residual))*1e6);
        assert(bound*1e6<=cfg.maxExtrapolationErrorUs, ...
            '%s C%02d: extrapolation uncertainty %.2f us exceeds limit.',s,prn,bound*1e6);
        [~,da,db]=intersect(epochKeys(ad.rows(:,1)),epochKeys(bd.rows(:,1)));
        good=bd.rows(db,4)>=3; da=da(good); db=db(good);
        n=numel(da); errors=sum(ad.rows(da,2)~=bd.rows(db,2));
        assert(n>=cfg.minimumMatchedBits && errors/max(n,1)<=cfg.maxBitErrorRate, ...
            '%s C%02d: air-bit match failed (%d/%d).',s,prn,errors,n);
        data=transfer(ad,ao,offset,scale,cfg.B.durationSeconds);
        overlay=transfer(ao,ao,offset,scale,cfg.B.durationSeconds);
        [covered,gaps]=coverage(data,cfg.B.durationSeconds);
        if ~covered && cfg.requireFullCoverage
            error(['%s C%02d: observed A air bits do not cover B=[0,%.3f) s. ' ...
                'Uncovered B intervals: %s. Start A earlier/end later or shorten B.'], ...
                s,prn,cfg.B.durationSeconds,mat2str(gaps,6));
        end
        if ~covered, warning('%s C%02d uncovered B intervals %s.',s,prn,mat2str(gaps,6)); end
        rec=struct('PRN',prn,'offsetAminusBSeconds',offset,'scaleAperB',scale, ...
            'clockDriftPpm',(scale-1)*1e6,'rmsResidualUs',sqrt(mean(residual.^2))*1e6, ...
            'maxResidualUs',max(abs(residual))*1e6,'extrapolationBoundUs',bound*1e6, ...
            'matchedBits',n,'bitErrors',errors,'fullCoverage',covered,'gaps',gaps,'method',method);
        if strcmp(s,'B1C')
            rec.secondaryBlockCoherence=bo.phaseCheck.coherence;
            rec.secondaryGlobalCoherence=bo.phaseCheck.globalCoherence;
            rec.secondaryGoodBlockFraction=bo.phaseCheck.goodBlockFraction;
            rec.dataPilotPhaseConsistency=bo.phaseCheck.dataPilotSquaredCoherence;
            rec.pilotPolarityTransitions=bo.phaseCheck.phasePolarityTransitions;
        end
        if isempty(records), records=rec; else, records(end+1)=rec; end %#ok<AGROW>
        allOutput{end+1}=struct('signal',s,'PRN',prn,'data',data,'overlay',overlay, ...
            'prefixData',bd.rows,'prefixOverlay',bo.rows, ...
            'week',ad.week,'step',ad.step,'overlayStep',ao.step,'fit',rec); %#ok<AGROW>
        fprintf('%s C%02d: A=%.9f + %.12f*B; residual %.3f us; air bits %d/%d errors; coverage %d.\n', ...
            s,prn,offset,scale,rec.maxResidualUs,errors,n,covered);
    end
    alignment.(s)=records;
end
common=intersect([alignment.B1I.PRN],[alignment.B1C.PRN]);
assert(~isempty(common),'No common B1I/B1C aids.');
for prn=common
    i=alignment.B1I([alignment.B1I.PRN]==prn); c=alignment.B1C([alignment.B1C.PRN]==prn);
    assert(abs(i.offsetAminusBSeconds-c.offsetAminusBSeconds)<50e-6, ...
        'C%02d: cross-signal A/B offsets disagree by >50 us.',prn);
end
for k=1:numel(allOutput)
    o=allOutput{k};
    for scope={'B_aligned','B_prefix'}
        out=fullfile(cfg.outputDir,scope{1},o.signal); if ~isfolder(out), mkdir(out); end
        base=sprintf('%s_PRN%02d',o.signal,o.PRN);
        if strcmp(scope{1},'B_aligned')
            data=o.data; overlay=o.overlay; kind='A_AIR_DATA_ALIGNED_TO_B';
        else, data=o.prefixData; overlay=o.prefixOverlay; kind='OBSERVED_B_PREFIX_AIR_BITS'; end
        writeRows(fullfile(out,[base '_data_bits.txt']),data,o,cfg,o.step,kind);
        if strcmp(o.signal,'B1I'), suffix='_NH.txt'; else, suffix='_secondary_code.txt'; end
        writeRows(fullfile(out,[base suffix]),overlay,o,cfg,o.overlayStep,'SYNCHRONIZED_OVERLAY');
    end
end
alignment.commonPRNs=common;
end

function [bd,bo]=matchB1IBits(o,ad,ao,cfg)
b=o.data; use=find(b(:,4)>=3 & b(:,6)<=cfg.B.strongSeconds);
assert(numel(use)>=cfg.minimumMatchedBits,'Not enough B air bits.');
observed=2*b(use,2)-1; reference=2*ad.rows(:,2)-1;
assert(numel(reference)>=numel(observed),'A reference is shorter than B prefix.');
scores=conv(reference,flipud(observed),'valid');
[peak,start]=max(abs(scores)); errors=(numel(observed)-peak)/2;
assert(errors/numel(observed)<=cfg.maxBitErrorRate,'B1I air-bit sequence not found in A.');
other=abs(scores); other(start)=-Inf;
assert(peak-max(other)>=6,'Ambiguous A/B air-bit match. Increase strongSeconds.');
if scores(start)<0, b(:,2)=1-b(:,2); b(:,3)=-b(:,3); end
first=start-use(1)+1; index=first+(0:size(b,1)-1)';
keep=index>=1 & index<=size(ad.rows,1); b=b(keep,:); index=index(keep);
b(:,1)=ad.rows(index,1); b(:,5)=ad.rows(index,5);
bd=ad; bd.rows=b;
bo=ao; bo.rows=o.overlay;
anchor=find(keep,1); anchorLoop=o.startIndices(anchor);
anchorSow=ad.rows(index(1),1);
bo.rows(:,1)=mod(anchorSow+((1:size(bo.rows,1))'-anchorLoop)*0.001,604800);
[~,ia,ib]=intersect(epochKeys(ao.rows(:,1)),epochKeys(bo.rows(:,1)));
assert(all(ao.rows(ia,5)==bo.rows(ib,5)),'A/B NH phase disagreement.');
end
function [bd,bo]=matchB1C(t,s,ad,ao,model,cfg)
rx=t.absoluteSample(:)/s.samplingFreq-t.remCodePhase(:)./t.codeFreq(:) ...
    -s.filterDelaySeconds-s.segmentStartSecond;
predicted=model.offsetAminusBSeconds+model.scaleAperB*rx;
index=interp1(ao.rows(:,4),(1:size(ao.rows,1))',predicted,'nearest',nan);
assert(all(isfinite(index)),'A secondary reference does not cover B prefix.');
assert(max(abs(ao.rows(index,4)-predicted))<50e-6,'B1I/B1C primary epoch mismatch.');
bo=ao; bo.rows=ao.rows(index,:); bo.rows(:,4)=rx;
bo.rows(:,6)=[rx(2:end);rx(end)+s.intTime];
use=rx>=3;
pilot=complex(t.Pilot_I_P(:),t.Pilot_Q_P(:));
z=complex(t.I_P(:),t.Q_P(:));
[bits,phaseCheck]=resolveB1CCarrierPolarity(z,pilot,bo.rows(:,3),rx,s.intTime,use);
assert(phaseCheck.valid, ...
    ['C%02d: secondary/phase check failed: block coherence %.3f, global %.3f, ' ...
     'good blocks %.3f, data/pilot phase consistency %.3f.'], ...
    t.PRN,phaseCheck.coherence,phaseCheck.globalCoherence, ...
    phaseCheck.goodBlockFraction,phaseCheck.dataPilotSquaredCoherence);
bo.phaseCheck=phaseCheck;
[~,aa,bb]=intersect(epochKeys(ad.rows(:,1)),epochKeys(bo.rows(:,1)));
good=rx(bb)>=3; aa=aa(good); bb=bb(good);
err=sum(bits(bb)~=ad.rows(aa,2));
if err>numel(bb)/2, bits=1-bits; err=numel(bb)-err; end
assert(numel(bb)>=cfg.minimumMatchedBits && err/numel(bb)<=cfg.maxBitErrorRate, ...
    'C%02d: B1C air-bit verification failed (%d/%d).',t.PRN,err,numel(bb));
bd=ad; bd.rows=[bo.rows(:,1) bits 1-2*bits rx bo.rows(:,5:6) zeros(numel(rx),1)];
fprintf('C%02d: B1C secondary block coherence %.3f (global %.3f); phase polarity handled; air-bit errors %d/%d.\n', ...
    t.PRN,phaseCheck.coherence,phaseCheck.globalCoherence,err,numel(bb));
end

function key=epochKeys(t)
key=int64(round(mod(t,604800)*1000));
end
function t=unwrapSow(t)
t=t(:); t=t+604800*cumsum([0;diff(t)<-302400]);
end
function [a,b,r,bound]=clockFit(ta,tb,cfg)
% Per-second medians reduce serially correlated DLL noise in drift fitting.
groups=floor(tb); g=unique(groups); x=zeros(size(g)); y=x;
for k=1:numel(g), ix=groups==g(k); x(k)=median(tb(ix)); y(k)=median(ta(ix)); end
if strcmpi(cfg.clockModel,'offset')
    b=1; a=median(ta-tb);
elseif strcmpi(cfg.clockModel,'affine')
    assert(numel(g)>=10,'Affine clock fit needs >=10 one-second windows.');
    p=polyfit(x-x(1),y-y(1),1); b=p(1); a=y(1)+p(2)-b*x(1);
    assert(abs(b-1)<100e-6,'Unreasonable A/B clock drift (>100 ppm).');
else, error('clockModel must be offset or affine.'); end
r=ta-a-b*tb;
% Conservative stability check from the two halves of the strong prefix.
mid=floor(numel(x)/2);
p1=polyfit(x(1:mid),y(1:mid),1); p2=polyfit(x(mid+1:end),y(mid+1:end),1);
horizon=max(0,cfg.B.durationSeconds-max(tb));
if b==1, drift=abs(polyfitSlope(x,y)-1); else, drift=abs(p1(1)-p2(1)); end
bound=max(abs(r))+drift*horizon;
end
function slope=polyfitSlope(x,y)
p=polyfit(x-x(1),y-y(1),1); slope=p(1);
end
function rows=transfer(aid,overlay,offset,scale,duration)
rows=aid.rows;
% Interpolate code epochs rather than adding nominal dt to sample times.
ot=unwrapSow(overlay.rows(:,1)); dt=unwrapSow(rows(:,1));
dt=dt+604800*round((ot(1)-dt(1))/604800);
ends=interp1(ot,overlay.rows(:,4),dt+aid.step,'linear','extrap');
rows(:,4)=(rows(:,4)-offset)/scale;
rows(:,6)=(ends-offset)/scale;
rows=rows(rows(:,6)>0 & rows(:,4)<duration,:);
end
function [ok,gaps]=coverage(rows,duration)
gaps=[]; cursor=0; tolerance=2e-6;
for k=1:size(rows,1)
    if rows(k,4)>cursor+tolerance, gaps(end+1,:)=[cursor,min(rows(k,4),duration)]; end %#ok<AGROW>
    cursor=max(cursor,rows(k,6));
end
if cursor<duration-tolerance, gaps(end+1,:)=[cursor duration]; end
ok=isempty(gaps);
end
function writeRows(file,rows,o,cfg,step,kind)
fid=fopen(file,'wt','n','UTF-8'); assert(fid>0,'Cannot write %s.',file);
closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'# %s PRN%02d %s; B-relative air bits/codes\n',o.signal,o.PRN,kind);
fprintf(fid,'# BDT_WEEK=%g; bit_duration_seconds=%.3f; PM1=1-2*BIT\n',o.week,step);
fprintf(fid,'# B_SOURCE=%s\n',char(java.io.File(cfg.B.file).getCanonicalPath()));
info=readLabsat3wInfo(cfg.B.file);
fprintf(fid,'# B_FILE_BYTES=%.0f; B_FILE_DATENUM=%.12f\n',info.fileBytes,info.modifiedDatenum);
fprintf(fid,'# B_FILE_START_SECONDS=%.9f; B_DURATION_SECONDS=%.9f\n',cfg.B.startSecond,cfg.B.durationSeconds);
fprintf(fid,'# A_RX=%.12f+%.15f*B_RX; residual_max_us=%.6f; extrapolation_bound_us=%.6f\n', ...
    o.fit.offsetAminusBSeconds,o.fit.scaleAperB,o.fit.maxResidualUs,o.fit.extrapolationBoundUs);
fprintf(fid,'# RX_SECONDS/END_SECONDS: B-relative receive interval [start,end); PHASE_INDEX is 1-based.\n');
fprintf(fid,'# First row may start below 0: this is the bit/chip spanning B start. No unknown bits are filled.\n');
fprintf(fid,'# FEC_OK=1: checked/re-encoded; 0: raw air-bit decision. No complete-frame output gate.\n');
fprintf(fid,'# For OVERLAY, FEC_OK=1 denotes synchronized deterministic code phase.\n');
if size(rows,2)<7, rows(:,7)=1; end
fprintf(fid,'BDT_SOW\tBIT\tPM1\tRX_SECONDS\tPHASE_INDEX\tEND_SECONDS\tFEC_OK\n');
fprintf(fid,'%.3f\t%d\t%d\t%.9f\t%d\t%.9f\t%d\n',rows(:,1:7).');
end

