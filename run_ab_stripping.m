function result = run_ab_stripping(cfg)
% A strong reference -> B strong prefix -> aligned aids for ALL of B.
% No B samples after strongSeconds are used for navigation alignment.
root=fileparts(mfilename('fullpath')); addpath(root,fullfile(root,'tools'));
if nargin<1, cfg=ab_stripping_settings(); end
assert(isfile(cfg.A.file) && isfile(cfg.B.file),'Set cfg.A.file and cfg.B.file.');
aInfo=readLabsat3wInfo(cfg.A.file); bInfo=readLabsat3wInfo(cfg.B.file);
if isempty(cfg.B.durationSeconds)
    cfg.B.durationSeconds=bInfo.durationSeconds-cfg.B.startSecond;
end
assert(cfg.A.startSecond>=0 && cfg.B.startSecond>=0,'Negative segment start.');
assert(cfg.A.processSeconds>0 && cfg.B.strongSeconds>=12 && ...
    cfg.B.strongSeconds<=cfg.B.durationSeconds,'Invalid A/B processing duration.');
for t=[cfg.A.processSeconds cfg.B.strongSeconds]
    assert(abs(t*100-round(t*100))<1e-6,'Processing lengths must be multiples of 10 ms.');
end
assert(cfg.A.startSecond+cfg.A.processSeconds+0.011<aInfo.durationSeconds, ...
    'A needs at least 11 ms beyond the requested tracking interval (code offset).');
assert(cfg.B.startSecond+cfg.B.durationSeconds<=bInfo.durationSeconds+1e-9,'B window exceeds recording.');
if strcmpi(char(java.io.File(cfg.A.file).getCanonicalPath()), ...
        char(java.io.File(cfg.B.file).getCanonicalPath()))
    fprintf('SAME-RECORDING DEMO: B=[%.3f,%.3f) file seconds. Tracking/alignment use the original strong samples.\n', ...
        cfg.B.startSecond,cfg.B.startSecond+cfg.B.durationSeconds);
end
if ~isfolder(cfg.cacheDir), mkdir(cfg.cacheDir); end
% Separate source caches prevent equal-length A/B files overwriting each other.
pa=convert_labsat3w(cfg.A.file,conversionLength(cfg.A.startSecond+cfg.A.processSeconds,aInfo), ...
    fullfile(cfg.cacheDir,'IQ_A'));
if isequal(aInfo,bInfo) && cfg.B.startSecond+cfg.B.strongSeconds+0.02<=pa.B1I.seconds
    pb=pa;
else
    pb=convert_labsat3w(cfg.B.file,conversionLength(cfg.B.startSecond+cfg.B.strongSeconds,bInfo), ...
        fullfile(cfg.cacheDir,'IQ_B'));
end
prefix=cfg.B; prefix.processSeconds=cfg.B.strongSeconds;
for signal={'B1I','B1C'}
    s=signal{1};
    result.A.(s)=receiverSegment(root,pa,s,cfg.A,cfg,'A');
end
% Report the ORIGINAL strong A before B processing / weak acquisition.
% Recomputed from prompts even when matching tracking caches are reused.
result.A.cn0=exportStrongReferenceCN0(result.A,cfg);
for signal={'B1I','B1C'}
    s=signal{1};
    result.B.(s)=receiverSegment(root,pb,s,prefix,cfg,'B_prefix');
end
result.alignment=alignAidsToB(result.A,result.B,cfg);
result.cfg=cfg;
% Metadata and tracking are intermediates, outside the data-only output tree.
save(fullfile(cfg.cacheDir,'aligned_run.mat'),'result','-v7.3');
if cfg.runAcquisition
    addpath(fullfile(root,'LS3W_Real_B1I_B1C'));
    acq=cfg.acquisition;
    acq.ls3wFile=cfg.B.file; acq.fileStartSecond=cfg.B.startSecond;
    acq.durationSeconds=cfg.B.durationSeconds;
    acq.aidDir=fullfile(cfg.outputDir,'B_aligned');
    acq.prnList=result.alignment.commonPRNs;
    acq.referenceRunFile=fullfile(cfg.cacheDir,'aligned_run.mat');
    result.acquisition=Joint_B1I_B1C_PhaseComp_Acquisition(acq);
    save(fullfile(cfg.cacheDir,'aligned_run.mat'),'result','-v7.3');
end
end
function n=conversionLength(last,info)
% Use 10 ms chunks; leave room for code phase and the receiver read buffer.
n=min(ceil((last+1)*100)/100,floor(info.durationSeconds*100)/100);
end
