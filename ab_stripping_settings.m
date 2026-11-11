function cfg = ab_stripping_settings()
% Code Creator: Taibin Liu    Code Creation Date: 2026.08
% A/B processing configuration. Times are seconds relative to each FILE.
% The shipped configuration is an explicit SAME-RECORDING functional demo.
% result = run_ab_stripping;
root=fileparts(mfilename('fullpath'));
cfg.A.file=fullfile(root,'data','File_048.LS3W');
cfg.A.startSecond=0;
cfg.A.processSeconds=80;
% 每次A处理后（包括复用缓存）打印并保存强信号载噪比统计。
cfg.A.cn0WarmupSeconds=5;         % 跳过初始环路收敛段；相对A处理起点。
cfg.A.cn0WindowSeconds=1;         % 每个统计窗口的长度。
cfg.A.cn0MinimumCoherence=0.5;    % 剔除载波相干性不足的窗口。
cfg.A.cn0ResultsDir=fullfile(root,'reports','A强信号载噪比');
cfg.B.file=cfg.A.file;            % Replace with your actual B recording.   cfg.A.file
cfg.B.startSecond=25;             % B=[25,80)秒；B相对40秒对应原文件65秒。
cfg.B.durationSeconds=55;         % [] = remaining duration of B.file.
cfg.B.strongSeconds=25;           % ONLY this prefix is acquired/tracked.
cfg.prns=7;                      % 同录制示例为C07；当前NH辅助模式支持6:58。
cfg.maxChannels=12;               % 实际进入跟踪的通道上限，避免全PRN均开跟踪。
% A/B原接收机的初始捕获范围；与后面的弱信号联合搜索网格独立。
cfg.receiverAcquisition.B1I.searchBandHz=5e3; % 单边范围，即[-100,+100] kHz。
cfg.receiverAcquisition.B1I.frequencyStepHz=500;
cfg.receiverAcquisition.B1I.noncoherentMs=20;
cfg.receiverAcquisition.B1C.searchBandHz=5e3;
cfg.receiverAcquisition.B1C.frequencyStepHz=50; % 10ms主码捕获需保留足够细的频率步长。
cfg.clockModel='offset';          % Shared clock. 'affine' estimates clock drift.
cfg.maxAlignmentResidualUs=20;
cfg.maxExtrapolationErrorUs=50;    % Reject uncertain long extrapolation.
cfg.maxBitErrorRate=0.02;
cfg.minimumMatchedBits=60;
cfg.requireFullCoverage=true;     % Do not silently fill unvalidated nav gaps.
cfg.reuseTracking=true;
cfg.outputDir=fullfile(root,'output','AB_File048_demo');
cfg.cacheDir=fullfile(root,'cache','ab_stripping');
cfg.runAcquisition=true;
cfg.acquisition.startSecond=40;   % Relative to B WINDOW start, not A/file.
cfg.acquisition.showAReference=true; % 捕获结束后仅在命令窗口对照A跟踪及约200ms短捕获，不保存参考结果。
% 仅在本次捕获读取的 B 窗口加噪；A 和 B 前 strongSeconds 秒保持原样。
% cn0DropDb 是 C/N0 的下降量(dB)，不是目标 C/N0，也不是幅度倍数。
cfg.acquisition.noise.enabled=false; % 公开默认不加噪；弱信号试验设true。
cfg.acquisition.noise.cn0DropDb=34; % 下降量：47 dB-Hz下降35 dB后约12 dB-Hz。
cfg.acquisition.noise.seed=2026; % 相同文件绝对样点和种子复现相同噪声。
cfg.acquisition.noise.calibrationMs=50; % 在捕获起点读取原始强信号，估计两支路噪声底。
cfg.acquisition.noise.calibrationMode='code-spectrum'; % 码频谱加权等效噪声底；'bandpower'保留旧的全带均值法对照。
cfg.acquisition.noise.originalNoisePsd=[]; % []自动估计；或手动填[B1I B1C]复基带双边PSD，原始解包单位^2/Hz。
cfg.acquisition.noise.referenceCN0dBHz=[]; % 可选：已测得的[B1I B1C]原始C/N0，仅记录理论下降后数值，不参与加噪。
cfg.acquisition.b1iCoherentMs=10; % 公开默认10ms；弱信号实验可改为700等。
cfg.acquisition.b1cCoherentMs=10; % 可与B1I不同，例如730ms。
% B1C alone/joint share receive windows and C background at corresponding
% phase/Doppler hypotheses; standalone still searches all ten primary slots.
cfg.acquisition.jointWeights=[1 1]; % [B1I B1C] COMPLEX-AMPLITUDE weights.
cfg.acquisition.jointWeightNormalization='l2'; % Scale invariant; [4 1] == [.8 .2].  none
% 跨B1I/B1C的相对载波相位补偿（不是B1C内部数据/导频的固定90度）。
cfg.acquisition.jointPhaseMode='continuous'; % 'continuous': 闭式连续最优；'grid': 离散相位。
cfg.acquisition.jointPhaseBins=360; % 仅grid生效：K个点覆盖[0,360)，间隔360/K度。
% 例：K=16 -> 22.5度；K=360 -> 1度；K=3600 -> 0.1度。
% For 80/20 percent POWER coefficients, use jointWeights=sqrt([.8 .2]).
% Any positive duration, including fractional ms; quantized only to samples.
cfg.acquisition.autoDopplerStep=true; % Refine frequency step to <=0.5/T.
cfg.acquisition.codeDriftCompensation=true; % A-derived code-epoch time warp.
cfg.acquisition.useGPU='auto';        % 'auto', true, or false (CPU fallback).
cfg.acquisition.dopplerRateHzPerSecond=0; % Optional known carrier Doppler rate.
cfg.acquisition.dopplerHz=-500:10:500; % B1C参考网格，自动细化只减小间隔。
% 当前10ms积分可用50Hz步长；长积分自动细化，宽频段耗时会随之增加。
cfg.acquisition.b1cBranch='both'; % 'both': data + pilot BOC11; or 'pilot'/'data'.
cfg.acquisition.b1cDataPilotWeights=sqrt([11 29]); % [data pilot11] amplitude weights; L2 normalized.
cfg.acquisition.b1cPilotPhaseDeg=90; % SIGNAL pilot11 phase relative to data; correlator rotates by -90 deg.
% 三个门限独立：单路是否通过，不影响联合判决。
cfg.acquisition.b1iMinRatio=8;    % 仅用于 B1I 单独捕获：B1IRatio >= 此值。
cfg.acquisition.b1cMinRatio=8;    % 仅用于 B1C 单独捕获：B1CRatio >= 此值。
cfg.acquisition.minJointScore=16; % 仅用于联合捕获：jointScore >= 此值。
cfg.acquisition.figureVisible='on';
cfg.acquisition.saveFigures=true;
cfg.acquisition.figureLayout='separate'; % 'separate':两张1x3；'six':一张2x3，上三维下侧视，每列一个底标。
cfg.acquisition.resultsDir=fullfile(root,'results','AB_acquisition');
% 独立的 C/N0 加噪线性验证：运行 validate_cn0_attenuation。
% 不运行长积分联合搜索；固定强信号码/载波参考，重新解扩原始/加噪样本。
cfg.cn0Validation.prn=7;
cfg.cn0Validation.startSecond=cfg.acquisition.startSecond; % B 相对时间
cfg.cn0Validation.durationSeconds=10;
cfg.cn0Validation.attenuationDb=0:1:15; % 各分支仅取预计C/N0 >= 下限+余量的点
cfg.cn0Validation.minimumCN0dBHz=35;
cfg.cn0Validation.marginDb=0.5;
cfg.cn0Validation.windowSeconds=1; % 用于估计结果的分窗离散程度
cfg.cn0Validation.b1cPromptMs=1; % 1或10；强信号矩估计默认统一1ms，保留10ms诊断
cfg.cn0Validation.seed=cfg.acquisition.noise.seed;
cfg.cn0Validation.dopplerHz=[]; % []从A强信号跟踪缓存取得；或[B1I B1C] Hz
cfg.cn0Validation.reusePrompts=true;
cfg.cn0Validation.figureVisible='on';
cfg.cn0Validation.resultsDir=fullfile(root,'reports','CN0加噪线性验证');
end
