function results = Joint_B1I_B1C_PhaseComp_Acquisition(cfg)
% 从 B 的可配置时刻读取实测 LS3W，使用 A/B 对齐后的 TXT 辅助捕获。
% 先运行根目录 run_ab_stripping；参数统一见 ab_stripping_settings.m。
% 不传参时使用该配置；也可传入 cfg.acquisition 及下列输入字段。
% B1I: 电文*NH；B1C: pilot/data/both，both按固定正交相位合成两路。
% 已知 NH 相位直接决定 B1C 的 1 ms 搜索槽，不要求弱 B1I 先捕获成功。
% B1I/B1C/联合完整三维搜索图保存在 results；output 仍然只有码/电文。
% 捕获窗口可选数字加噪见 cfg.noise；以统一配置为准，强前缀不加噪。
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'tools'));
base=ab_stripping_settings(); defaults=base.acquisition;
defaults.ls3wFile=base.B.file;
defaults.fileStartSecond=base.B.startSecond;
defaults.durationSeconds=base.B.durationSeconds;
defaults.aidDir=fullfile(base.outputDir,'B_aligned');
defaults.prnList=base.prns;
defaults.referenceRunFile=fullfile(base.cacheDir,'aligned_run.mat');
if nargin<1, cfg=struct(); end
if isfield(cfg,'minBranchRatio')
    if ~isfield(cfg,'b1iMinRatio'), cfg.b1iMinRatio=cfg.minBranchRatio; end
    if ~isfield(cfg,'b1cMinRatio'), cfg.b1cMinRatio=cfg.minBranchRatio; end
end
% Explicit legacy calls using only coherentMs keep their original duration
% and unnormalized statistic rather than inheriting new branch defaults.
if isfield(cfg,'coherentMs') && ~isfield(cfg,'b1iCoherentMs') && ~isfield(cfg,'b1cCoherentMs')
    cfg.b1iCoherentMs=cfg.coherentMs; cfg.b1cCoherentMs=cfg.coherentMs;
    if ~isfield(cfg,'jointWeightNormalization'), cfg.jointWeightNormalization='none'; end
    if ~isfield(cfg,'minJointScore') && strcmpi(cfg.jointWeightNormalization,'none')
        cfg.minJointScore=32;
    end
end
for f=fieldnames(defaults).'
    if ~isfield(cfg,f{1}), cfg.(f{1})=defaults.(f{1}); end
end
if isempty(cfg.durationSeconds)
    info=readLabsat3wInfo(cfg.ls3wFile);
    cfg.durationSeconds=info.durationSeconds-cfg.fileStartSecond;
end
if isempty(cfg.prnList)
    files=dir(fullfile(cfg.aidDir,'B1I','*_NH.txt')); available=[];
    for k=1:numel(files)
        tok=regexp(files(k).name,'PRN(\d+)','tokens','once');
        prn=str2double(tok{1});
        if isfile(fullfile(cfg.aidDir,'B1C',sprintf('B1C_PRN%02d_secondary_code.txt',prn)))
            available(end+1)=prn; %#ok<AGROW>
        end
    end
    cfg.prnList=available;
end
assert(~isempty(cfg.prnList),'No common aligned aids; run run_ab_stripping first.');
% Snapshot the alignment metadata before a potentially long weak search.
% A reference captures run afterwards and are printed only, never saved.
referenceRun=[];
if cfg.showAReference && isfile(cfg.referenceRunFile)
    try
        saved=load(cfg.referenceRunFile,'result'); referenceRun=saved.result;
    catch err
        warning('AReference:Unavailable','Cannot load A reference run: %s',err.message);
    end
end
results=aidedJointAcquisition(cfg);
if cfg.showAReference
    if isempty(referenceRun)
        warning('AReference:Unavailable','A reference unavailable: run run_ab_stripping first (%s).',cfg.referenceRunFile);
    else
        try
            printAcquisitionAReference(referenceRun,cfg,results);
        catch err
            warning('AReference:Unavailable','B acquisition completed; A console reference unavailable: %s',err.message);
        end
    end
end
end
