function result = validate_cn0_attenuation(options)
%VALIDATE_CN0_ATTENUATION Measure, rather than impose, C/N0 after digital AWGN.
% result = validate_cn0_attenuation; % configuration: ab_stripping_settings
% options overrides cfg.cn0Validation fields only. No acquisition/tracking
% settings are changed. Existing A tracking and B_aligned epochs are needed.
% Complex prompts default to 1 ms for this strong-signal moment test;
% B1C can also use 10 ms. This does not change receiver/acquisition periods.
% Carrier phase need
% not be recovered for M2M4; within-prompt Doppler/code mismatch is fixed.
root=fileparts(mfilename('fullpath')); old=path; restore=onCleanup(@()path(old));
addpath(root,fullfile(root,'tools'),fullfile(root,'mBDS_B1I','include'), ...
    fullfile(root,'mBDS-3_B1C','include'));
base=ab_stripping_settings(); v=base.cn0Validation;
if nargin<1, options=struct(); end
assert(isstruct(options) && isscalar(options),'options must be a scalar struct.');
for f=fieldnames(options).'
    assert(isfield(v,f{1}),'Unknown validation option: %s',f{1});
    v.(f{1})=options.(f{1});
end
validateOptions(v,base);
work=fullfile(root,'cache','cn0_noise_validation');
if ~isfolder(work), mkdir(work); end
if ~isfolder(v.resultsDir), mkdir(v.resultsDir); end
c=base.acquisition; c.ls3wFile=base.B.file; c.fileStartSecond=base.B.startSecond;
c.startSecond=v.startSecond; c.noise.enabled=true; c.noise.cn0DropDb=10;
c.noise.seed=v.seed; c.noise.referenceCN0dBHz=[];
model=prepareAcquisitionNoise(c); % SAME automatic/manual calibration as acquisition
aidDir=fullfile(base.outputDir,'B_aligned');
fileI=fullfile(aidDir,'B1I',sprintf('B1I_PRN%02d_NH.txt',v.prn));
fileC=fullfile(aidDir,'B1C',sprintf('B1C_PRN%02d_secondary_code.txt',v.prn));
aidI=readStrippingAid(fileI); aidC=readStrippingAid(fileC);
checkAid(aidI,base); checkAid(aidC,base);
[fd,reference]=referenceDoppler(base,v,aidI,aidC);
key=struct('schema',2,'source',readLabsat3wInfo(base.B.file), ...
    'fileStartSecond',base.B.startSecond,'prn',v.prn,'start',v.startSecond, ...
    'seconds',v.durationSeconds,'dopplerHz',fd,'noiseModel',model, ...
    'aidIdentity',[dir(fileI) dir(fileC)]);
cacheFile=fullfile(work,sprintf('PRN%02d_prompts.mat',v.prn));
reuse=false;
if v.reusePrompts && isfile(cacheFile)
    saved=load(cacheFile,'key'); reuse=isequaln(saved.key,key);
end
if reuse
    saved=load(cacheFile,'prompts'); prompts=saved.prompts;
    fprintf('Using matching raw-domain prompt cache.\n');
else
    prompts=extractPrompts(base,v,model,aidI,aidC,fd);
    save(cacheFile,'key','prompts','reference','-v7.3');
end
% Every operation up to a fixed prompt is linear. The noise prompts below
% come from REAL raw-domain AWGN through the SAME reader/FIR/decimation/code
% correlators. Scaling these cached noise prompts is exactly equivalent to
% reprocessing the same seeded raw noise at every D, without repeated FIR.
% An independent direct raw read at D=3 verifies this equivalence.
if v.b1cPromptMs==1
    z0={prompts.I0,prompts.D10,prompts.P10};
    zn={prompts.IN,prompts.D1N,prompts.P1N};
    times={prompts.timeI,prompts.timeI,prompts.timeI};
else
    z0={prompts.I0,prompts.D0,prompts.P0};
    zn={prompts.IN,prompts.DN,prompts.PN};
    times={prompts.timeI,prompts.timeC,prompts.timeC};
end
period=[.001 v.b1cPromptMs/1000 v.b1cPromptMs/1000];
names={'B1I','B1C data','B1C pilot11','B1C total (data + pilot11)'};
baseline=zeros(1,4);
for s=1:3, baseline(s)=momentCN0(z0{s},period(s)); end
baseline(4)=powerSumDb(baseline(2),baseline(3));
assert(all(isfinite(baseline)),'Invalid original C/N0 estimate; inspect reference alignment.');
fprintf('Measured ORIGINAL C/N0 [I data pilot11 total]=%s dB-Hz.\n',mat2str(baseline,6));
drop=sort(unique(double(v.attenuationDb(:))));
drop=unique([0;drop]);
% The total is a sum of separately measured branches. Require BOTH B1C
% components above the requested floor instead of relying on a weak data
% estimate just because the total/pilot remains strong.
maximum=[baseline(1:3),min(baseline(2:3))]-v.minimumCN0dBHz-v.marginDb;
rows=struct([]); blockRows=struct([]);
for id=1:numel(drop)
    d=drop(id); scale=sqrt(expm1(log(10)*d/10)/expm1(log(10)));
    estimate=nan(1,4); block=nan(floor(v.durationSeconds/v.windowSeconds),4);
    for s=1:3
        if d>maximum(s)+1e-9, continue; end
        z=z0{s}+scale*zn{s}; estimate(s)=momentCN0(z,period(s));
        for j=1:size(block,1)
            lo=v.startSecond+(j-1)*v.windowSeconds;
            ix=times{s}>=lo & times{s}<lo+v.windowSeconds;
            block(j,s)=momentCN0(z(ix),period(s));
        end
    end
    if d<=maximum(4)+1e-9
        estimate(4)=powerSumDb(estimate(2),estimate(3));
        block(:,4)=powerSumDb(block(:,2),block(:,3));
    end
    for s=1:4
        if d>maximum(s)+1e-9, continue; end
        prediction=baseline(s)-d;
        eligible=isfinite(estimate(s)) && estimate(s)>=v.minimumCN0dBHz;
        row=struct('PRN',v.prn,'Branch',names{s},'ConfiguredDropDb',d, ...
            'PromptMilliseconds',period(min(s,3))*1000, ...
            'ValidationStartBSecond',v.startSecond,'ValidationDurationSeconds',v.durationSeconds, ...
            'CalibrationMode',model.calibrationMode, ...
            'OriginalCN0dBHz',baseline(s),'TheoreticalCN0dBHz',prediction, ...
            'MeasuredCN0dBHz',estimate(s),'MeasuredDropDb',baseline(s)-estimate(s), ...
            'ErrorVsTheoryDb',estimate(s)-prediction, ...
            'WindowStdDb',std(block(isfinite(block(:,s)),s),0), ...
            'ValidWindows',sum(isfinite(block(:,s))), ...
            'MeasuredAtLeastMinimum',eligible,'Seed',v.seed);
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
        for j=1:size(block,1)
            r=struct('Branch',names{s},'ConfiguredDropDb',d, ...
                'WindowStartBSecond',v.startSecond+(j-1)*v.windowSeconds, ...
                'MeasuredCN0dBHz',block(j,s),'Valid',isfinite(block(j,s)));
            if isempty(blockRows), blockRows=r; else, blockRows(end+1)=r; end %#ok<AGROW>
        end
    end
end
assert(~isempty(rows),'Original signal is below the requested validation floor.');
tableData=struct2table(rows); fits=fitBranches(rows,names);
legacyRows=legacyComparison(rows,z0,zn,period,baseline,model,names,v.minimumCN0dBHz);
if isempty(legacyRows), legacyFits=struct([]); else, legacyFits=fitBranches(legacyRows,names); end
result=struct('options',v,'baselineCN0dBHz',baseline,'rows',rows,'fits',fits, ...
    'legacyRows',legacyRows,'legacyFits',legacyFits, ...
    'noiseModel',model,'carrierReference',reference,'promptCounts', ...
    cellfun(@numel,z0), ...
    'directRawPromptRelativeError',prompts.directRawPromptRelativeError, ...
    'method','Raw-LS3W reference-aided prompts; M2M4; no per-D parameter fitting; B1C linear power sum');
writetable(tableData,fullfile(v.resultsDir,'cn0_attenuation_measurements.csv'),'Encoding','UTF-8');
writetable(struct2table(fits),fullfile(v.resultsDir,'cn0_attenuation_fit.csv'),'Encoding','UTF-8');
writetable(struct2table(blockRows),fullfile(v.resultsDir,'cn0_attenuation_windows.csv'),'Encoding','UTF-8');
writetable(struct2table(legacyRows),fullfile(v.resultsDir,'cn0_attenuation_legacy_bandpower.csv'),'Encoding','UTF-8');
save(fullfile(work,'validation_result.mat'),'result','-v7.3');
fig=drawFigure(result); result.figure=fig;
writeReport(result,v.resultsDir);
disp(tableData(:,{'Branch','ConfiguredDropDb','OriginalCN0dBHz','MeasuredCN0dBHz','ErrorVsTheoryDb'}));
disp(struct2table(fits));
fprintf('Validation finished: %s\n',v.resultsDir);
end

function rows=legacyComparison(rows,z0,zn,period,baseline,model,names,minimum)
if ~all(isfinite(model.bandpowerNoisePsd)), rows=rows([]); return; end
r=sqrt(model.bandpowerNoisePsd./model.originalNoisePsd);
for k=1:numel(rows)
    d=rows(k).ConfiguredDropDb; scale=sqrt(expm1(log(10)*d/10)/9);
    db=zeros(1,3);
    for s=1:3
        band=1+(s>1);
        db(s)=momentCN0(z0{s}+scale*r(band)*zn{s},period(s));
    end
    db(4)=powerSumDb(db(2),db(3));
    s=find(strcmp(rows(k).Branch,names),1);
    rows(k).MeasuredCN0dBHz=db(s);
    rows(k).MeasuredDropDb=baseline(s)-db(s);
    rows(k).ErrorVsTheoryDb=db(s)-(baseline(s)-d);
    rows(k).WindowStdDb=NaN; rows(k).ValidWindows=0;
    rows(k).CalibrationMode='bandpower';
    rows(k).MeasuredAtLeastMinimum=isfinite(db(s)) && db(s)>=minimum;
end
end

function p=extractPrompts(base,v,model,aidI,aidC,fd)
fs=model.rawFs/4;
begin=v.startSecond; finish=begin+v.durationSeconds;
rxI=selectedEpochs(aidI,begin,finish); rxC=selectedEpochs(aidC,begin,finish);
codeI=double(generateB12Icode(v.prn)); codeI=codeI(:);
s.codeLength=10230;
codeD=generateDataBOC11(s,v.prn); codeD=double(codeD(1:6:end)); codeD=codeD(:);
codeP=generatePilotBOC11(s,v.prn); codeP=double(codeP(1:6:end)); codeP=codeP(:);
sumI=complex(zeros(numel(rxI),2)); sumD=complex(zeros(numel(rxC),2)); sumP=sumD;
sumD1=sumI; sumP1=sumI;
countI=zeros(numel(rxI),1); countC=zeros(numel(rxC),1);
nTotal=round(v.durationSeconds*fs); first=round((base.B.startSecond+begin)*fs);
chunk=round(.1*fs); timer=tic; nextUpdate=0; directError=NaN;
for offset=0:chunk:nTotal-1
    n=min(chunk,nTotal-offset); absolute=(first+offset)/fs;
    [aI,aC]=readLabsatWindow(base.B.file,absolute,n/fs*1000);
    [bI,bC]=readLabsatWindow(base.B.file,absolute,n/fs*1000,model);
    xI=[aI,bI-aI]; xC=[aC,bC-aC];
    t=(offset+(0:n-1)')/fs+begin;
    % Fixed carrier reference; no carrier/noise/signal fit on degraded data.
    mixI=exp(-2j*pi*mod((t-begin)*fd(1),1));
    mixC=exp(-2j*pi*mod((t-begin)*fd(2),1));
    [epochI,chipI]=codePosition(rxI,t,2046);
    [epochC,chipC]=codePosition(rxC,t,20460);
    refI=mixI.*codeI(chipI); refD=mixC.*codeD(chipC); refP=mixC.*codeP(chipC);
    for col=1:2
        sumI(:,col)=sumI(:,col)+groupSum(epochI,xI(:,col).*refI,numel(rxI));
        sumD(:,col)=sumD(:,col)+groupSum(epochC,xC(:,col).*refD,numel(rxC));
        sumP(:,col)=sumP(:,col)+groupSum(epochC,xC(:,col).*refP,numel(rxC));
        sumD1(:,col)=sumD1(:,col)+groupSum(epochI,xC(:,col).*refD,numel(rxI));
        sumP1(:,col)=sumP1(:,col)+groupSum(epochI,xC(:,col).*refP,numel(rxI));
    end
    countI=countI+accumarray(epochI,ones(n,1),[numel(rxI) 1]);
    countC=countC+accumarray(epochC,ones(n,1),[numel(rxC) 1]);
    if offset==0
        direct=model; direct.requestedDropDb=3; direct.appliedDropDb=3;
        direct.addedNoisePsd=model.originalNoisePsd*expm1(log(10)*.3);
        [dI,dC]=readLabsatWindow(base.B.file,absolute,n/fs*1000,direct);
        scale=sqrt(expm1(log(10)*.3)/expm1(log(10)));
        actual=[sum(dI.*refI),sum(dC.*refD),sum(dC.*refP)];
        predicted=[sum((aI+scale*xI(:,2)).*refI), ...
            sum((aC+scale*xC(:,2)).*refD),sum((aC+scale*xC(:,2)).*refP)];
        directError=max(abs(actual-predicted)./max(1,abs(actual)));
        assert(directError<1e-8,'Cached linear prompts disagree with direct raw-domain noise.');
        fprintf('Direct D=3 dB raw-domain equivalence: relative prompt error %.3g.\n',directError);
    end
    if toc(timer)>=nextUpdate || offset+n==nTotal
        fprintf('Raw-domain C/N0 prompts: %.1f/%.1f s processed, %.1f s elapsed.\n', ...
            (offset+n)/fs,v.durationSeconds,toc(timer));
        nextUpdate=toc(timer)+15;
    end
end
fullI=rxI(1:end-1)>=begin & rxI(2:end)<=finish;
fullC=rxC(1:end-1)>=begin & rxC(2:end)<=finish;
ix=find(fullI); cx=find(fullC);
assert(all(countI(ix)>.999*.001*fs) && all(countC(cx)>.999*.01*fs),'Incomplete prompt accumulation.');
p=struct('I0',sumI(ix,1)./countI(ix),'IN',sumI(ix,2)./countI(ix), ...
    'D0',sumD(cx,1)./countC(cx),'DN',sumD(cx,2)./countC(cx), ...
    'P0',sumP(cx,1)./countC(cx),'PN',sumP(cx,2)./countC(cx), ...
    'D10',sumD1(ix,1)./countI(ix),'D1N',sumD1(ix,2)./countI(ix), ...
    'P10',sumP1(ix,1)./countI(ix),'P1N',sumP1(ix,2)./countI(ix), ...
    'timeI',(rxI(ix)+rxI(ix+1))/2,'timeC',(rxC(cx)+rxC(cx+1))/2, ...
    'directRawPromptRelativeError',directError);
end

function rx=selectedEpochs(aid,a,b)
rx=aid.rows(:,4);
lo=max(1,find(rx<=a,1,'last')-1); hi=find(rx>=b,1,'first')+1;
assert(~isempty(lo) && ~isempty(hi) && hi<=numel(rx),'Insufficient aligned epoch coverage.');
rx=rx(lo:hi);
end
function [epoch,chip]=codePosition(rx,t,lengthCode)
u=interp1(rx,(1:numel(rx))',t,'linear');
assert(all(isfinite(u)),'Code phase beyond reference epochs.');
epoch=floor(u); frac=u-epoch;
chip=min(lengthCode,floor(frac*lengthCode)+1);
end
function z=groupSum(group,x,n)
z=accumarray(group,real(x),[n 1])+1j*accumarray(group,imag(x),[n 1]);
end
function value=momentCN0(z,T)
value=NaN;
if numel(z)<30 || any(~isfinite(z)), return; end
power=abs(z).^2; m2=mean(power); rad=2*m2^2-mean(power.^2);
if rad<=0 || m2<=0, return; end
c=sqrt(rad); noise=m2-c;
if noise<=0, return; end
value=10*log10(c/noise/T);
end
function db=powerSumDb(a,b)
db=10*log10(10.^(a/10)+10.^(b/10));
end

function [fd,ref]=referenceDoppler(base,v,aidI,aidC)
if ~isempty(v.dopplerHz)
    fd=double(v.dopplerHz(:).'); ref=struct('method','manual fixed Doppler','Hz',fd); return;
end
fd=zeros(1,2); aids={aidI,aidC};
for k=1:2
    signals={'B1I','B1C'};
    file=fullfile(base.cacheDir,'A',[signals{k} '_tracking.mat']);
    assert(isfile(file),'Run run_ab_stripping first, or supply cn0Validation.dopplerHz.');
    q=load(file,'trackResults','settings');
    ch=find([q.trackResults.PRN]==v.prn,1); assert(~isempty(ch),'A tracking lacks selected PRN.');
    tok=regexp(aids{k}.header,'A_RX=([-+\d.eE]+)\+([-+\d.eE]+)\*B_RX','tokens','once');
    assert(~isempty(tok),'Missing A/B receive-time mapping.');
    offset=str2double(tok{1}); scale=str2double(tok{2}); t=q.trackResults(ch);
    rx=t.absoluteSample(:)/q.settings.samplingFreq-q.settings.filterDelaySeconds-q.settings.segmentStartSecond;
    ia=rx>=offset+scale*v.startSecond & rx<offset+scale*(v.startSecond+v.durationSeconds);
    assert(sum(ia)>100,'A tracking does not cover validation window.');
    f=double(t.carrFreq(ia))-q.settings.IF;
    fd(k)=median(f)*scale;
    assert(max(abs(f-median(f)))<25, ...
        'Carrier dynamics exceed fixed-reference validation; choose a shorter window.');
end
ref=struct('method','median A strong tracking Doppler, fixed before adding noise','Hz',fd);
fprintf('Fixed strong-reference Doppler [I C]=[%.6f %.6f] Hz.\n',fd);
end
function checkAid(aid,base)
info=readLabsat3wInfo(base.B.file);
source=regexp(aid.header,'(?m)^# B_SOURCE=([^\r\n]*)','tokens','once');
assert(~isempty(source) && strcmpi(strtrim(source{1}),char(java.io.File(base.B.file).getCanonicalPath())), ...
    'Aid belongs to another B recording.');
keys={'B_FILE_BYTES','B_FILE_DATENUM','B_FILE_START_SECONDS','B_DURATION_SECONDS'};
expected=[info.fileBytes,info.modifiedDatenum,base.B.startSecond,base.B.durationSeconds];
tol=[0 1e-9 1/info.SMP 1e-6];
for k=1:4
    tok=regexp(aid.header,[keys{k} '=([-+\d.eE]+)'],'tokens','once');
    assert(~isempty(tok) && abs(str2double(tok{1})-expected(k))<=tol(k),'Stale/mismatched B aid: %s',keys{k});
end
end
function validateOptions(v,base)
assert(isscalar(v.prn) && ismember(v.prn,6:58),'Select a D1 PRN 6:58.');
assert(isscalar(v.startSecond) && isfinite(v.startSecond) && v.startSecond>=0,'Invalid validation start.');
assert(isscalar(v.durationSeconds) && isfinite(v.durationSeconds) && v.durationSeconds>=2 && ...
    v.startSecond+v.durationSeconds+0.02<base.B.durationSeconds,'Invalid validation duration.');
assert(isvector(v.attenuationDb) && ~isempty(v.attenuationDb) && ...
    all(isfinite(v.attenuationDb)) && all(v.attenuationDb>=0) && all(v.attenuationDb<=60),'Invalid attenuation grid.');
assert(isscalar(v.minimumCN0dBHz) && isfinite(v.minimumCN0dBHz) && v.minimumCN0dBHz>=35,'Validation minimum must be >=35 dB-Hz.');
assert(isscalar(v.marginDb) && isfinite(v.marginDb) && v.marginDb>=0,'Invalid margin.');
assert(isscalar(v.windowSeconds) && v.windowSeconds>=.5 && v.windowSeconds<=v.durationSeconds/2,'Invalid statistical window.');
assert(isscalar(v.b1cPromptMs) && ismember(v.b1cPromptMs,[1 10]),'b1cPromptMs must be 1 or 10.');
assert(isempty(v.dopplerHz) || (isreal(v.dopplerHz) && numel(v.dopplerHz)==2 && all(isfinite(v.dopplerHz))),'dopplerHz must be [] or [B1I B1C].');
end
function fits=fitBranches(rows,names)
fits=struct([]);
for s=1:4
    % Keep every point selected BEFORE seeing noisy estimates. Do not make
    % the fitted trend look better by censoring deviations below the floor.
    a=rows(strcmp({rows.Branch},names{s}));
    assert(all(isfinite([a.MeasuredCN0dBHz])),'Invalid M2M4 estimate for %s.',names{s});
    assert(numel(a)>=3,'Fewer than three strong-range points for %s.',names{s});
    x=[a.ConfiguredDropDb]; y=[a.MeasuredCN0dBHz]; fit=polyfit(x,y,1);
    err=y-([a.OriginalCN0dBHz]-x); pred=polyval(fit,x);
    r=struct('Branch',names{s},'PointCount',numel(a),'MaximumTestedDropDb',max(x), ...
        'BaselineCN0dBHz',a(1).OriginalCN0dBHz,'SlopeDbPerDb',fit(1), ...
        'InterceptdBHz',fit(2),'R2',1-sum((y-pred).^2)/sum((y-mean(y)).^2), ...
        'MaxAbsErrorVsTheoryDb',max(abs(err)),'RMSErrorVsTheoryDb',sqrt(mean(err.^2)));
    if isempty(fits), fits=r; else, fits(end+1)=r; end %#ok<AGROW>
end
end

function fig=drawFigure(result)
v=result.options; fits=result.fits;
colors=[.16 .35 .55;.72 .40 .19;.13 .52 .48;.36 .36 .58];
fig=figure('Color','w','Position',[100 70 1100 800],'Visible',v.figureVisible,'InvertHardcopy','off');
set(fig,'PaperUnits','inches','PaperPosition',[0 0 7.2 5.4],'PaperSize',[7.2 5.4]);
tl=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for s=1:4
    ax=nexttile(tl); hold(ax,'on');
    a=result.rows(strcmp({result.rows.Branch},fits(s).Branch));
    x=[a.ConfiguredDropDb]; y=[a.MeasuredCN0dBHz]; theory=[a.TheoreticalCN0dBHz];
    yline(ax,fits(s).BaselineCN0dBHz,':','Color',[.57 .60 .64],'LineWidth',1.0,'DisplayName','Original baseline');
    plot(ax,x,theory,'--','Color',[.20 .22 .25],'LineWidth',1.2,'DisplayName','Baseline - set drop');
    errorbar(ax,x,y,[a.WindowStdDb],'-o','Color',colors(s,:),'MarkerFaceColor','w', ...
        'LineWidth',1.25,'MarkerSize',4,'CapSize',3,'DisplayName','Measured (M2M4)');
    ax.FontName='Arial'; ax.FontSize=9; ax.LineWidth=.7; ax.TickDir='out'; ax.Box='off';
    ax.Color='w'; ax.YGrid='on'; ax.GridAlpha=.13; ax.XGrid='off';
    xlim(ax,[0 max(x)+.15]);
    spread=[a.WindowStdDb]; spread(~isfinite(spread))=0;
    % Include the COMPLETE theoretical line and error bars, even when a
    % biased 10 ms estimate sits above theory. Never clip contrary evidence.
    low=min([y-spread,theory])-.4;
    high=max([y+spread,theory,fits(s).BaselineCN0dBHz])+.6;
    ylim(ax,[low high]);
    xlabel(ax,'Configured C/N_0 reduction (dB)'); ylabel(ax,'Estimated C/N_0 (dB-Hz)');
    title(ax,sprintf('%c  %s',char('a'+s-1),fits(s).Branch),'FontWeight','bold','FontSize',10);
    text(ax,.97,.80,sprintf('Baseline %.2f dB-Hz\nSlope %.3f; max error %.2f dB', ...
        fits(s).BaselineCN0dBHz,fits(s).SlopeDbPerDb,fits(s).MaxAbsErrorVsTheoryDb), ...
        'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
        'FontSize',8,'Color',[.25 .28 .32]);
    if s==1, legend(ax,'Location','southwest','Box','off','FontSize',7.8); end
end
title(tl,sprintf('Raw-sample noise validation | PRN %02d | B %.0f-%.0f s', ...
    v.prn,v.startSecond,v.startSecond+v.durationSeconds),'FontName','Arial','FontSize',12,'FontWeight','normal');
subtitle(tl,sprintf('Fixed reference; 1 ms B1I / %g ms B1C prompts; error bars: %.0f s-window SD', ...
    v.b1cPromptMs,v.windowSeconds),'FontName','Arial','FontSize',8,'Color',[.35 .38 .42]);
drawnow;
file=fullfile(v.resultsDir,'cn0_attenuation_validation');
print(fig,[file '.png'],'-dpng','-r600'); print(fig,[file '.pdf'],'-dpdf','-painters');
print(fig,[file '.svg'],'-dsvg','-painters'); savefig(fig,[file '.fig']);
end
function writeReport(result,folder)
fid=fopen(fullfile(folder,'验证结果说明.md'),'wt','n','UTF-8'); assert(fid>0,'Cannot write report.');
closer=onCleanup(@()fclose(fid));
v=result.options;
fprintf(fid,'# 原始样本加噪后的载噪比线性验证\n\n');
fprintf(fid,'PRN %02d；B [%.3f, %.3f) 秒；同一录制片段、同一噪声种子 %u。固定使用强信号码边界和载波频率参考，所有衰减档位采用相同相关器；未重新优化弱信号码/载波参数。\n\n', ...
    v.prn,v.startSecond,v.startSecond+v.durationSeconds,v.seed);
fprintf(fid,'| 分支 | 原始估计(dB-Hz) | 测试下降范围(dB) | 拟合斜率(dB/dB) | R² | 相对理论线最大偏差(dB) |\n| --- | --- | --- | --- | --- | --- |\n');
for r=result.fits
    fprintf(fid,'| %s | %.3f | 0–%g | %.4f | %.6f | %.3f |\n',r.Branch,r.BaselineCN0dBHz,r.MaximumTestedDropDb,r.SlopeDbPerDb,r.R2,r.MaxAbsErrorVsTheoryDb);
end
fprintf(fid,'\n理论线为原始实测基准减去配置下降量，理想斜率是 -1；实测点来自重新解扩后复提示相关值的 M2M4 估计，未用该理论线生成实测点。\n\n');
fprintf(fid,'## 原噪声标定方法的对照\n\n');
fprintf(fid,'本轮发现原录制B1C背景谱不平坦，旧的全带平均功率法不足以代表码相关器见到的噪声。现默认使用独立的Welch谱估计，按BPSK(2)/BOC(1,1)码频谱加权标定等效噪声PSD；没有用M2M4目标值反推加噪幅度。旧方法对照保存在 `cn0_attenuation_legacy_bandpower.csv`。\n\n');
fprintf(fid,'| 分支 | 旧全带均值法斜率 | 当前码频谱加权斜率 |\n| --- | --- | --- |\n');
for k=1:numel(result.legacyFits)
    fprintf(fid,'| %s | %.4f | %.4f |\n',result.fits(k).Branch,result.legacyFits(k).SlopeDbPerDb,result.fits(k).SlopeDbPerDb);
end
fprintf(fid,'\n全带均值等效PSD [I C]=%s，码频谱加权等效PSD=%s，单位为解包数值平方/Hz。不能再仅凭全带功率是否增加D dB判定单星C/N0是否下降D dB。\n\n', ...
    mat2str(result.noiseModel.bandpowerNoisePsd,6),mat2str(result.noiseModel.codeWeightedNoisePsd,6));
fprintf(fid,'M2=mean(|z|²)，M4=mean(|z|⁴)，C=sqrt(2M2²-M4)，N=M2-C，C/N0=10log10(C/(N*T))。B1I 的 T=1 ms，B1C 数据和导频11各为%g ms；合成值=10log10(10^(数据/10)+10^(导频/10))，不补偿未匹配的BOC(6,1)。\n\n',v.b1cPromptMs);
fprintf(fid,'只在基准减去衰减不低于 %.1f dB-Hz 的区间取点（下限 %.1f 加 %.1f dB 余量）；不按结果好坏删点，若实测低于下限会在 CSV 标明。合成值要求数据、导频两者均在验证范围。误差条是相同录制内 %.1f 秒分窗估计的标准差，不是独立噪声重复试验的置信区间。拟合与R²仅描述本轮曲线，不等同于独立重复试验统计检验。\n\n', ...
    v.minimumCN0dBHz+v.marginDb,v.minimumCN0dBHz,v.marginDb,v.windowSeconds);
fprintf(fid,'噪声是在实际 LS3W 解包后、各路原 FIR 之前加入；原始与加噪样本均重新解扩。为避免每档重复处理同一录制，分别缓存原始相关值和经相同滤波/解扩得到的新增噪声相关值，再按原始噪声幅度缩放；所有这些处理均是固定线性运算。另用3 dB档直接从原始样本加噪重算，复相关值相对差 %.3g，检查了与完整原始样本路径的一致性。\n\n', ...
    result.directRawPromptRelativeError);
fprintf(fid,'该测试验证数字加噪及其载噪比变化，不评价弱信号重新捕获/跟踪的能力。原有信号幅度波动、码/载波参考误差、M2M4有限样本偏差和自动噪声底标定误差均会造成偏离；未将原始载噪比假设为45 dB-Hz，也未把带内总功率变化当作单星载噪比。\n\n');
if v.b1cPromptMs==1
    fprintf(fid,'另附 [B1C 10 ms相关周期对照](B1C_10ms对照/验证结果说明.md)。1 ms和10 ms估计口径必须区分；主图的误差不能直接外推到10 ms矩估计或原接收机经过环路、平滑后的界面数值。10 ms对照使用同一原始片段和同一噪声，只改变统计前的相关周期，未重新调节噪声标定来强行对齐曲线。\n\n');
end
fprintf(fid,'程序：根目录 `validate_cn0_attenuation.m`；参数：`ab_stripping_settings.m` 的 `cfg.cn0Validation`。直接运行 `result=validate_cn0_attenuation;`。图、四张数据CSV和本说明位于本目录；相关值缓存位于 `cache/cn0_noise_validation`，output 不增加文件。\n\n');
fprintf(fid,'方法定义参考 [GNSS-SDR 载噪比估计](https://gnss-sdr.org/docs/sp-blocks/tracking/#carrier-to-noise-density-ratio)。\n');
end
