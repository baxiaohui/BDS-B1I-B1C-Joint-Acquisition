function result = compare_b1i_b1c(phase, signal)
%COMPARE_B1I_B1C Fresh acquisition/tracking and matched-time C/N0 analysis.
% compare_b1i_b1c runs acquisition, both trackers, then statistics.
% Intermediate MAT/JSON files stay in cache/cno_power_comparison.
% Phase calls allow independent B1I/B1C tracking in separate MATLAB sessions:
% compare_b1i_b1c('acquire'); compare_b1i_b1c('track','B1I');
% compare_b1i_b1c('track','B1C'); compare_b1i_b1c('analyze');
if nargin<1, phase='run'; end
if nargin<2, signal=''; end
root=fileparts(mfilename('fullpath'));
work=fullfile(root,'cache','cno_power_comparison');
if ~isfolder(work), mkdir(work); end
basePath=path; cleanup=onCleanup(@()path(basePath));
addpath(root,fullfile(root,'tools'));
result=[];
if strcmp(phase,'run')
    compare_b1i_b1c('acquire');
    compare_b1i_b1c('track','B1I');
    compare_b1i_b1c('track','B1C');
    result=compare_b1i_b1c('analyze'); return;
end
if strcmp(phase,'acquire')
    config=struct('sourceFile',fullfile(root,'data','File_048.LS3W'), ...
        'msToProcess',45000,'warmupSeconds',5,'windowSeconds',1, ...
        'minimumCoherence',0.5,'sourceInfo',[], ...
        'runStarted',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
    converted=convert_labsat3w(config.sourceFile,46);
    acquired=cell(1,2);
    for k=1:2
        if k==1, s='B1I'; receiver='mBDS_B1I'; else, s='B1C'; receiver='mBDS-3_B1C'; end
        path(basePath); addpath(root,fullfile(root,'tools'),fullfile(root,receiver), ...
            fullfile(root,receiver,'include'),fullfile(root,receiver,'Common'));
        clear initSettings acquisitionCPU acquisitionGPU;
        settings=initSettings(); p=converted.(s);
        settings.fileName=p.fileName; settings.inputFormat='IQ_BIN';
        settings.dataType=p.dataType; settings.fileType=p.fileType;
        settings.samplingFreq=p.samplingFreq; settings.IF=p.IF;
        settings.originalFileName=p.source.fileName;
        settings.sourceInfo=p.source; settings.filterDelaySeconds=p.filterDelaySeconds;
        settings.msToProcess=config.msToProcess;
        settings.correlatorType=0; settings.trkMode=1;
        settings.plotTracking=0; settings.showWaitbar=false;
        settings.saveNavMessages=false;
        config.sourceInfo=p.source;
        fid=fopen(settings.fileName,'rb'); assert(fid>0,'Cannot open converted IQ.');
        closeInput=onCleanup(@()fclose(fid));
        fprintf('Fresh %s acquisition started.\n',s); tic;
        if exist('gpuDeviceCount','file') && gpuDeviceCount>0
            acqResults=acquisitionGPU(fid,settings); settings.gpuACQflag=1;
        else
            acqResults=acquisitionCPU(fid,settings); settings.gpuACQflag=0;
        end
        fprintf('Fresh %s acquisition elapsed %.1f s.\n',s,toc);
        clear closeInput;
        acquired{k}=find(acqResults.carrFreq~=0);
        save(fullfile(work,[s '_acquisition.mat']),'settings','acqResults');
    end
    config.commonPRNs=intersect(acquired{1},acquired{2});
    config.B1IAcquiredPRNs=acquired{1}; config.B1CAcquiredPRNs=acquired{2};
    assert(~isempty(config.commonPRNs),'No common acquired PRN.');
    save(fullfile(work,'config.mat'),'config');
    fprintf('Jointly acquired PRNs: %s\n',mat2str(config.commonPRNs));
    result=config; return;
end
load(fullfile(work,'config.mat'),'config');
if strcmp(phase,'track')
    assert(ismember(string(signal),["B1I","B1C"]),'Choose B1I or B1C.');
    if strcmp(signal,'B1I'), receiver='mBDS_B1I'; else, receiver='mBDS-3_B1C'; end
    addpath(fullfile(root,receiver),fullfile(root,receiver,'include'),fullfile(root,receiver,'Common'));
    clear preRun trkChannelsParallel corrMatlabParallelBPSK corrMatlabParallelB1C;
    load(fullfile(work,[signal '_acquisition.mat']),'settings','acqResults');
    assert(isequal(settings.sourceInfo,readLabsat3wInfo(config.sourceFile)), ...
        'Source recording changed since acquisition.');
    acqResults.carrFreq(setdiff(1:numel(acqResults.carrFreq),config.commonPRNs))=0;
    settings.numberOfChannels=numel(config.commonPRNs);
    channel=preRun(acqResults,settings);
    showChannelStatus(channel,settings);
    fid=fopen(settings.fileName,'rb'); assert(fid>0,'Cannot open converted IQ.');
    closeInput=onCleanup(@()closeIfOpen(fid)); tic;
    [trackResults,channel]=trkChannelsParallel(fid,channel,settings);
    elapsedSeconds=toc;
    assert(all([trackResults.status]=='T'),'Tracking ended early.');
    save(fullfile(work,[signal '_tracking.mat']), ...
        'settings','acqResults','channel','trackResults','elapsedSeconds');
    fprintf('Fresh %s tracking complete: %.1f s wall time.\n',signal,elapsedSeconds);
    result=trackResults; return;
end
assert(strcmp(phase,'analyze'),'Unknown phase.');
i=load(fullfile(work,'B1I_tracking.mat'));
c=load(fullfile(work,'B1C_tracking.mat'));
rows=struct('PRN',{},'startSeconds',{},'endSeconds',{},'B1ICNo',{}, ...
    'B1CDataCNo',{},'B1CPilotCNo',{},'B1ICoherence',{},'B1CDataCoherence',{}, ...
    'B1CPilotCoherence',{},'B1ICount',{},'B1CCount',{},'valid',{},'reason',{});
summary=struct('PRN',{},'validWindows',{},'totalWindows',{}, ...
    'B1ICNo',{},'B1CDataCNo',{},'B1CPilotCNo',{},'B1CTotalCNo',{}, ...
    'deltaCMinusI',{},'ratioCOverI',{},'ratioIOverC',{}, ...
    'deltaP10',{},'deltaP90',{},'pilotDataGap',{}, ...
    'B1IValidFrames',{},'B1CValidFrames',{},'status',{});
for prn=config.commonPRNs
    ti=i.trackResults(find([i.trackResults.PRN]==prn,1));
    tc=c.trackResults(find([c.trackResults.PRN]==prn,1));
    timeI=ti.absoluteSample/i.settings.samplingFreq-i.settings.filterDelaySeconds;
    timeC=tc.absoluteSample/c.settings.samplingFreq-c.settings.filterDelaySeconds;
    begin=numel(rows)+1;
    for a=config.warmupSeconds:config.windowSeconds:config.msToProcess/1000-config.windowSeconds
        b=a+config.windowSeconds;
        ix=find(timeI>=a & timeI+i.settings.intTime<=b);
        cx=find(timeC>=a & timeC+c.settings.intTime<=b);
        [gi,li,ri]=estimate(ti.I_P(ix),ti.Q_P(ix),i.settings.intTime);
        [gd,ld,rd]=estimate(tc.I_P(cx),tc.Q_P(cx),c.settings.intTime);
        [gp,lp,rp]=estimate(tc.Pilot_I_P(cx),tc.Pilot_Q_P(cx),c.settings.intTime);
        reason='OK';
        if numel(ix)<0.9*config.windowSeconds/i.settings.intTime || ...
                numel(cx)<0.9*config.windowSeconds/c.settings.intTime
            reason='INSUFFICIENT_SAMPLES';
        elseif ~isempty(ri), reason=['B1I_' ri];
        elseif ~isempty(rd), reason=['B1C_DATA_' rd];
        elseif ~isempty(rp), reason=['B1C_PILOT_' rp];
        elseif min([li ld lp])<config.minimumCoherence, reason='LOW_CARRIER_COHERENCE'; end
        rows(end+1)=struct('PRN',prn,'startSeconds',a,'endSeconds',b, ...
            'B1ICNo',gi,'B1CDataCNo',gd,'B1CPilotCNo',gp,'B1ICoherence',li, ...
            'B1CDataCoherence',ld,'B1CPilotCoherence',lp, ...
            'B1ICount',numel(ix),'B1CCount',numel(cx), ...
            'valid',strcmp(reason,'OK'),'reason',reason); %#ok<AGROW>
    end
    windowRows=rows(begin:end); eligible=windowRows([windowRows.valid]);
    n=numel(eligible); total=numel(windowRows);
    [framesI,framesC]=checkFrames(root,ti,c.trackResults,prn,c.settings);
    if n>0
        gi=linearMean([eligible.B1ICNo]); gd=linearMean([eligible.B1CDataCNo]);
        gp=linearMean([eligible.B1CPilotCNo]);
        gt=10*log10(10^(gd/10)+10^(gp/10));
        delta=gt-gi; ratio=10^(delta/10);
        combined=10*log10(10.^([eligible.B1CDataCNo]/10)+10.^([eligible.B1CPilotCNo]/10));
        differences=combined-[eligible.B1ICNo];
        p10=percentileLinear(differences,0.1); p90=percentileLinear(differences,0.9);
        gap=gp-gd;
        if n<0.8*total, status='LOW_VALID_COVERAGE';
        elseif framesI==0 || framesC==0, status='NO_FRAME_CONFIRMATION';
        else, status='FRAME_CONFIRMED'; end
    else
        gi=nan; gd=nan; gp=nan; gt=nan; delta=nan; ratio=nan;
        p10=nan; p90=nan; gap=nan; status='NO_RELIABLE_PAIRED_WINDOWS';
    end
    summary(end+1)=struct('PRN',prn,'validWindows',n,'totalWindows',total, ...
        'B1ICNo',gi,'B1CDataCNo',gd,'B1CPilotCNo',gp,'B1CTotalCNo',gt, ...
        'deltaCMinusI',delta,'ratioCOverI',ratio,'ratioIOverC',1/ratio, ...
        'deltaP10',p10,'deltaP90',p90,'pilotDataGap',gap, ...
        'B1IValidFrames',framesI,'B1CValidFrames',framesC,'status',status); %#ok<AGROW>
end
result=struct('config',config,'summary',summary,'windows',rows, ...
    'B1ITrackingWallSeconds',i.elapsedSeconds,'B1CTrackingWallSeconds',c.elapsedSeconds, ...
    'assumption','Equal input noise PSD; uncalibrated filter/quantization/tracking losses. B1C effective C/N0 is the linear sum of measured data and BOC(1,1) pilot C/N0; no missing-band extrapolation.', ...
    'method','Fresh tracking; paired receiver-time windows; M2M4 population moments; 10log10 of linear C/N0 means.');
save(fullfile(work,'comparison.mat'),'result');
fid=fopen(fullfile(work,'comparison.json'),'wt','n','UTF-8');
assert(fid>0,'Cannot write comparison JSON.'); closeJson=onCleanup(@()fclose(fid));
fprintf(fid,'%s',jsonencode(result)); clear closeJson;
disp(struct2table(summary));
end

function [ni,nc]=checkFrames(root,trackI,tracksC,prn,settingsC)
old=path; cleanup=onCleanup(@()path(old));
addpath(fullfile(root,'mBDS_B1I','include'),fullfile(root,'mBDS_B1I','Common'));
clear eph_structure_init ephemeris;
[~,~,~,messages]=NAVdecoding(trackI.I_P,prn); ni=sum([messages.valid]);
addpath(fullfile(root,'mBDS-3_B1C','include'),fullfile(root,'mBDS-3_B1C','Common'));
clear eph_structure_init ephemeris;
channel=find([tracksC.PRN]==prn,1);
[~,~,~,messages]=BCNAV1decoding(tracksC,channel,settingsC); nc=sum([messages.valid]);
end

function [cn0,coherence,reason]=estimate(I,Q,T)
cn0=nan; coherence=nan; reason='';
z=complex(double(I(:)),double(Q(:)));
if numel(z)<10 || any(~isfinite(z)), reason='INVALID_PROMPTS'; return; end
power=abs(z).^2; m2=mean(power); v=var(power,1);
if m2<=0, reason='ZERO_POWER'; return; end
coherence=abs(mean(z.^2))/m2;
radicand=m2^2-v;
if radicand<=0, reason='INVALID_MOMENT_RADICAND'; return; end
signalPower=sqrt(radicand); noisePower=m2-signalPower;
if noisePower<=0, reason='NONPOSITIVE_NOISE_POWER'; return; end
cn0=10*log10(signalPower/noisePower/T);
if ~isfinite(cn0), reason='NONFINITE_CNO'; end
end

function value=linearMean(db)
value=10*log10(mean(10.^(db/10)));
end

function value=percentileLinear(values,p)
values=sort(values); position=1+(numel(values)-1)*p;
lo=floor(position); hi=ceil(position);
value=values(lo)+(position-lo)*(values(hi)-values(lo));
end

function closeIfOpen(fid)
if ~isempty(fopen(fid)), fclose(fid); end
end
