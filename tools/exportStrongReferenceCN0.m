function report = exportStrongReferenceCN0(A,cfg)
%EXPORTSTRONGREFERENCECN0 Print/save measured C/N0 after strong-A processing.
% Uses un-smoothed ORIGINAL tracking prompts at the actual receiver periods,
% never the weak-acquisition surface or its artificial noise configuration.
root=fileparts(fileparts(mfilename('fullpath')));
options=struct('cn0WarmupSeconds',5,'cn0WindowSeconds',1, ...
    'cn0MinimumCoherence',.5,'cn0ResultsDir',fullfile(root,'reports','A强信号载噪比'));
for f=fieldnames(options).'
    if isfield(cfg.A,f{1}), options.(f{1})=cfg.A.(f{1}); end
end
warmup=options.cn0WarmupSeconds; width=options.cn0WindowSeconds;
coherenceMinimum=options.cn0MinimumCoherence;
assert(isscalar(warmup) && isfinite(warmup) && warmup>=0 && warmup<cfg.A.processSeconds, ...
    'A.cn0WarmupSeconds must lie in [0,A.processSeconds).');
assert(isscalar(width) && isfinite(width) && width>0 && width<=cfg.A.processSeconds-warmup, ...
    'A.cn0WindowSeconds must fit after the warmup interval.');
assert(isscalar(coherenceMinimum) && isfinite(coherenceMinimum) && ...
    coherenceMinimum>=0 && coherenceMinimum<=1,'A.cn0MinimumCoherence must be in [0,1].');
source=readLabsat3wInfo(cfg.A.file);
signals={'B1I','B1C'}; tracks=cell(1,2); settings=cell(1,2); acquisition=cell(1,2);
for k=1:2
    q=load(A.(signals{k}).trackingFile,'key','settings','trackResults','acqResults');
    assert(isequaln(q.key.source,source) && q.key.start==cfg.A.startSecond && ...
        q.key.seconds==cfg.A.processSeconds && isequal(q.key.prns,cfg.prns), ...
        'A C/N0 tracking cache does not match source/start/duration/PRNs: %s',signals{k});
    settings{k}=q.settings; tracks{k}=q.trackResults([q.trackResults.status]=='T');
    acquisition{k}=q.acqResults;
end
prns=union([tracks{1}.PRN],[tracks{2}.PRN]);
assert(~isempty(prns),'No tracked A satellites available for C/N0 statistics.');
count=floor((cfg.A.processSeconds-warmup)/width+1e-10);
windows=struct([]); summary=struct([]);
for prn=prns
    ti=tracks{1}([tracks{1}.PRN]==prn); tc=tracks{2}([tracks{2}.PRN]==prn);
    ai=acquisitionResult(acquisition{1},settings{1},prn);
    ac=acquisitionResult(acquisition{2},settings{2},prn);
    first=numel(windows)+1;
    for j=1:count
        a=warmup+(j-1)*width; b=a+width;
        [gi,li,ri,ni]=estimateWindow(ti,settings{1},a,b,'I_P','Q_P',coherenceMinimum);
        [gd,ld,rd,nc]=estimateWindow(tc,settings{2},a,b,'I_P','Q_P',coherenceMinimum);
        [gp,lp,rp]=estimateWindow(tc,settings{2},a,b,'Pilot_I_P','Pilot_Q_P',coherenceMinimum);
        validI=isempty(ri); validC=isempty(rd) && isempty(rp);
        gt=NaN; if validC, gt=powerSum(gd,gp); end
        row=struct('PRN',prn,'StartASecond',a,'EndASecond',b, ...
            'StartFileSecond',cfg.A.startSecond+a,'EndFileSecond',cfg.A.startSecond+b, ...
            'B1ICN0dBHz',gi,'B1CDataCN0dBHz',gd,'B1CPilot11CN0dBHz',gp, ...
            'B1CTotalCN0dBHz',gt,'B1ICoherence',li,'B1CDataCoherence',ld, ...
            'B1CPilotCoherence',lp,'B1IPromptCount',ni,'B1CPromptCount',nc, ...
            'B1IValid',validI,'B1CValid',validC,'PairedValid',validI && validC, ...
            'B1IStatus',statusText(ri),'B1CDataStatus',statusText(rd),'B1CPilotStatus',statusText(rp));
        if isempty(windows), windows=row; else, windows(end+1)=row; end %#ok<AGROW>
    end
    w=windows(first:end); wi=w([w.B1IValid]); wc=w([w.B1CValid]);
    gi=linearMean([wi.B1ICN0dBHz]); gd=linearMean([wc.B1CDataCN0dBHz]);
    gp=linearMean([wc.B1CPilot11CN0dBHz]); gt=powerSum(gd,gp);
    status='OK';
    if isempty(wi) || isempty(wc), status='NO_VALID_WINDOWS_FOR_ONE_OR_BOTH_SIGNALS';
    elseif min(numel(wi),numel(wc))<.8*count, status='LOW_VALID_COVERAGE'; end
    row=struct('PRN',prn,'AFileStartSecond',cfg.A.startSecond, ...
        'AProcessSeconds',cfg.A.processSeconds,'StatsStartASecond',warmup, ...
        'StatsEndASecond',warmup+count*width,'WindowSeconds',width, ...
        'B1IPromptMs',settings{1}.intTime*1000,'B1CPromptMs',settings{2}.intTime*1000, ...
        'B1IAcqDetected',ai.detected,'B1IAcqCodePhaseSample1',ai.codeSample1, ...
        'B1IAcqCodeDelayMs',ai.codeDelayMs,'B1IAcqDopplerHz',ai.dopplerHz, ...
        'B1IAcqCarrierHz',ai.carrierHz,'B1IAcqPeakMetric',ai.peakMetric, ...
        'B1IAcqZeroDopplerSentinel',ai.zeroDopplerSentinel, ...
        'B1CAcqDetected',ac.detected,'B1CAcqCodePhaseSample1',ac.codeSample1, ...
        'B1CAcqCodeDelayMs',ac.codeDelayMs,'B1CAcqDopplerHz',ac.dopplerHz, ...
        'B1CAcqCarrierHz',ac.carrierHz,'B1CAcqPeakMetric',ac.peakMetric, ...
        'B1CAcqZeroDopplerSentinel',ac.zeroDopplerSentinel, ...
        'B1ICN0dBHz',gi,'B1CDataCN0dBHz',gd,'B1CPilot11CN0dBHz',gp, ...
        'B1CTotalCN0dBHz',gt,'B1IValidWindows',numel(wi),'B1CValidWindows',numel(wc), ...
        'PairedValidWindows',sum([w.PairedValid]),'TotalWindows',count,'Status',status);
    if isempty(summary), summary=row; else, summary(end+1)=row; end %#ok<AGROW>
end
out=options.cn0ResultsDir; if ~isfolder(out), mkdir(out); end
[~,stem]=fileparts(cfg.A.file);
prefix=fullfile(out,sprintf('%s_A_%gs_%gs_CN0',stem,cfg.A.startSecond,cfg.A.processSeconds));
report=struct('summary',struct2table(summary),'windows',struct2table(windows), ...
    'summaryFile',[prefix '_summary.csv'],'windowsFile',[prefix '_windows.csv'], ...
    'methodFile',[prefix '_说明.md'],'options',options,'source',source, ...
    'method','Unsmoothed A prompts; M2M4 population moments; linear-domain mean; pilot BOC11 only');
report.acquisition=report.summary(:,{'PRN','AFileStartSecond', ...
    'B1IAcqDetected','B1IAcqCodePhaseSample1','B1IAcqCodeDelayMs','B1IAcqDopplerHz','B1IAcqPeakMetric', ...
    'B1CAcqDetected','B1CAcqCodePhaseSample1','B1CAcqCodeDelayMs','B1CAcqDopplerHz','B1CAcqPeakMetric'});
writetable(report.summary,report.summaryFile,'Encoding','UTF-8');
writetable(report.windows,report.windowsFile,'Encoding','UTF-8');
fprintf('\nA STRONG-REFERENCE C/N0 [dB-Hz] | A file %.3f s + %.3f s | statistics A [%.3f,%.3f) s\n', ...
    cfg.A.startSecond,cfg.A.processSeconds,warmup,warmup+count*width);
fprintf('Original prompts: B1I %g ms; B1C %g ms. M2M4; means in LINEAR domain.\n', ...
    settings{1}.intTime*1000,settings{2}.intTime*1000);
fprintf('A acquisition at A start: code delay=(stored sample index-1)/Fs; Doppler=carrier-IF [Hz].\n');
for k=1:height(report.acquisition)
    for s={'B1I','B1C'}
        signal=s{1};
        fprintf('  PRN%02d %s: codePhase=%.0f (1-based samples), codeDelay=%.9f ms, Doppler=%.3f Hz\n', ...
            report.acquisition.PRN(k),signal, ...
            report.acquisition.([signal 'AcqCodePhaseSample1'])(k), ...
            report.acquisition.([signal 'AcqCodeDelayMs'])(k), ...
            report.acquisition.([signal 'AcqDopplerHz'])(k));
    end
end
disp(report.summary(:,{'PRN','B1ICN0dBHz','B1CDataCN0dBHz','B1CPilot11CN0dBHz', ...
    'B1CTotalCN0dBHz','B1IValidWindows','B1CValidWindows','TotalWindows','Status'}));
fprintf('A C/N0 summary: %s\nA per-window C/N0: %s\n\n',report.summaryFile,report.windowsFile);
writeMethod(report,cfg);
end
function a=acquisitionResult(result,s,prn)
% CPU/GPU acquisition stores a 1-based complex-sample position. The
% tracking code uses skipNumberOfSamples+codePhase-1 for its first read.
a=struct('detected',false,'codeSample1',NaN,'codeDelayMs',NaN, ...
    'dopplerHz',NaN,'carrierHz',NaN,'peakMetric',NaN,'zeroDopplerSentinel',false);
if prn>numel(result.carrFreq) || ~isfinite(result.carrFreq(prn)) || result.carrFreq(prn)==0
    return;
end
codeSample=double(result.codePhase(prn));
assert(isfinite(codeSample) && codeSample>=1,'Invalid A acquisition code phase for PRN%02d.',prn);
a.detected=true; a.codeSample1=codeSample;
a.codeDelayMs=(codeSample-1)/s.samplingFreq*1000;
a.carrierHz=double(result.carrFreq(prn)); a.dopplerHz=a.carrierHz-s.IF;
% Original CPU/GPU fine acquisition reserves 1 Hz for a zero-frequency
% detection at zero IF (otherwise carrFreq==0 would mean not acquired).
a.zeroDopplerSentinel=s.IF==0 && a.carrierHz==1;
if a.zeroDopplerSentinel, a.dopplerHz=0; end
if isfield(result,'peakMetric') && numel(result.peakMetric)>=prn
    a.peakMetric=double(result.peakMetric(prn));
end
end
function [cn0,coherence,reason,n]=estimateWindow(t,s,a,b,iField,qField,minimum)
cn0=NaN; coherence=NaN; n=0; reason='NO_TRACKED_CHANNEL';
if isempty(t), return; end
time=double(t.absoluteSample(:))/s.samplingFreq-s.filterDelaySeconds-s.segmentStartSecond;
ix=find(time>=a & time+s.intTime<=b+1e-10); n=numel(ix);
if n<max(10,.9*(b-a)/s.intTime), reason='INSUFFICIENT_PROMPTS'; return; end
z=complex(double(t.(iField)(ix)),double(t.(qField)(ix))); z=z(:);
if any(~isfinite(z)), reason='NONFINITE_PROMPTS'; return; end
p=abs(z).^2; m2=mean(p);
if m2<=0, reason='ZERO_POWER'; return; end
coherence=abs(mean(z.^2))/m2;
rad=m2^2-var(p,1); % population moments: equivalent to 2*M2^2-M4
if rad<=0, reason='INVALID_MOMENT_RADICAND'; return; end
signal=sqrt(rad); noise=m2-signal;
if noise<=0, reason='NONPOSITIVE_NOISE'; return; end
cn0=10*log10(signal/noise/s.intTime);
if ~isfinite(cn0), reason='NONFINITE_CN0'; cn0=NaN; return; end
if coherence<minimum, reason='LOW_CARRIER_COHERENCE'; return; end
reason='';
end
function value=linearMean(db)
if isempty(db), value=NaN; else, value=10*log10(mean(10.^(db/10))); end
end
function value=powerSum(a,b)
value=10*log10(10.^(a/10)+10.^(b/10));
end
function s=statusText(reason)
if isempty(reason), s='OK'; else, s=reason; end
end
function writeMethod(report,cfg)
fid=fopen(report.methodFile,'wt','n','UTF-8'); assert(fid>0,'Cannot write A C/N0 description.');
closer=onCleanup(@()fclose(fid)); o=report.options;
fprintf(fid,'# A强信号载噪比输出\n\n');
fprintf(fid,'A文件：`%s`；从文件第 %.3f 秒起处理 %.3f 秒。跳过前 %.3f 秒，以 %.3f 秒为窗统计。\n\n', ...
    cfg.A.file,cfg.A.startSecond,cfg.A.processSeconds,o.cn0WarmupSeconds,o.cn0WindowSeconds);
fprintf(fid,'每次 `run_ab_stripping` 完成两路A处理后，自动打印并保存本表；复用跟踪缓存也会重新输出。汇总可通过 `result.A.cn0.summary` 读取，逐窗结果为 `result.A.cn0.windows`。\n\n');
fprintf(fid,'同一汇总CSV还保存A在处理起点处的B1I/B1C原始捕获结果，可用 `result.A.cn0.acquisition` 查看。AcqCodePhaseSample1是原捕获器保存的1基复采样点索引；AcqCodeDelayMs=(索引-1)/该路采样率*1000，多普勒AcqDopplerHz=捕获载波频率-IF。原始载波值另存AcqCarrierHz；AcqPeakMetric是原捕获器的主次峰比，与后续辅助联合捕获的背景归一化指标不同。未捕获的支路标为false及NaN。\n\n');
fprintf(fid,'码滞后相对A处理起点的转换后IQ输入窗口，B1I按1 ms主码周期、B1C按10 ms主码周期解释；这里保留捕获器坐标，不扣除转换FIR群延迟。AFileStartSecond给出捕获输入的文件起点。它不是整个A时段的平均码相位，也不是B后续时刻的码相位。零IF下原捕获器可能将零频检测存为1 Hz哨兵，报告恢复为0 Hz，并以AcqZeroDopplerSentinel标记；原存储值保留在AcqCarrierHz。\n\n');
fprintf(fid,'使用A未加噪的、未平滑的复提示相关值重新计算M2M4载噪比，B1I/B1C保持各自原接收机的相关周期（见表中PromptMs列），不等同于先前统一1 ms的独立加噪验证口径，也不直接平均界面平滑值。\n\n');
fprintf(fid,'M2=mean(|z|²)，M4=mean(|z|⁴)，S=sqrt(2*M2²-M4)，N=M2-S，C/N0=10*log10(S/(N*T))。汇总先在线性域求平均，再转为dB-Hz。C/N0由矩估计得到，表示当前接收机处理链的估计值。\n\n');
fprintf(fid,'B1C合成等效值=10*log10(10^(数据/10)+10^(导频11/10))。两路分量使用相同有效窗口汇总；导频为当前窄带BOC(1,1)分量，不外推未匹配的BOC(6,1)能量。B1I使用自身有效窗口，逐窗表保留PairedValid标识，便于选择跨信号共同有效时段。\n\n');
fprintf(fid,'提示相关值不足90%%、矩估计无效或载波相干性abs(mean(z²))/mean(|z|²)低于 %.3f 的窗口不计入均值；逐窗表保留状态，不静默填补。时间列StartASecond/EndASecond相对A处理窗口，StartFileSecond/EndFileSecond相对原始文件。\n\n',o.cn0MinimumCoherence);
fprintf(fid,'统计表保存在reports，output仍仅保存电文、NH和二级码TXT。修改A.cn0WarmupSeconds、A.cn0WindowSeconds可调整统计时段与分辨率。\n');
end
