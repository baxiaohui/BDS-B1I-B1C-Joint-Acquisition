北斗B1I与B1C强信号辅助联合捕获 使用说明

1. 功能与流程
LS3W读取与转换 → A强信号接收机处理 → B强前缀处理 → 电文、NH码、
二级码对齐 → 导出覆盖B的辅助符号 → B1I/B1C单路与加权相干联合捕获。
A应包含B的信号时间范围并留首尾裕量；未知电文不猜测补齐。

2. 运行准备
MATLAB R2024a；需要Signal Processing Toolbox、DSP System Toolbox、
Communications Toolbox。GPU可选，另需Parallel Computing Toolbox。
当前主流程可用CPU运行，不必编译MEX。
准备LS3W和同名INI，放入data目录，或在配置中填写自己的完整路径。
当前支持58 Msps、QUA=1、CHN=3、SFT=6，读取物理通道A。
源码不附录制或历史缓存，没有实际数据不能运行信号处理。

3. 怎么运行
在MATLAB中切换至本文件所在目录，执行：
  cfg = ab_stripping_settings();
  check_public_setup(cfg);       % 检查环境及数据，不做捕获跟踪
  result = run_public(cfg);      % 完整流程
也可直接执行：result = run_public;
不要对两套接收机同时addpath(genpath(...))，主流程会管理重名函数。

4. 在哪里改参数
统一修改ab_stripping_settings.m：
  A.file / startSecond / processSeconds：A路径、文件内起点、处理时长。
  B.file / startSecond / durationSeconds：B路径、文件内起点、窗口长度。
  B.strongSeconds：仅处理B的强前缀时长，当前至少12秒。
  prns：实际存在的同星PRN，NH辅助模式支持6至58。
  receiverAcquisition：强信号接收机的多普勒范围和步长。
  acquisition.startSecond：相对B窗口起点的捕获时刻。
  acquisition.b1iCoherentMs / b1cCoherentMs：两路独立积分毫秒数。
  acquisition.dopplerHz：B1C参考多普勒网格；B1I按载频比映射。
  acquisition.jointWeights：复幅度权重；jointWeightNormalization='l2'归一化。
  acquisition.jointPhaseMode：continuous闭式最优或grid离散相位网格。
  acquisition.b1cBranch：data、pilot或both；导频匹配BOC(1,1)部分。
  acquisition.b1iMinRatio / b1cMinRatio / minJointScore：独立检测门限。
  acquisition.figureLayout：separate两张三图；six一张2×3组合图。

默认示例：File_048、C07；A处理0至80秒；B取同录制25至80秒，
长度55秒，跟踪B前25秒，从B第40秒（原文件约65秒）捕获。
默认不加噪，两路积分10ms，强信号范围±5kHz，联合网格-500:10:500 Hz。
其他录制必须修改路径、PRN、多普勒及时间范围。A处理尾部之后还需
至少11ms以上真实录制。共用时钟用clockModel='offset'，独立时钟可用'affine'。

加噪实验：acquisition.noise.enabled=true；cn0DropDb是下降量而非目标。
例如原始47dB-Hz、cn0DropDb=35，理论上降至约12dB-Hz。
只给捕获使用的B窗口加噪；seed控制随机噪声，实际下降精度取决于标定。
弱信号可将两路积分改为700ms、730ms等，长积分会自动细化频率网格。

5. 主要代码
run_public.m：完整入口；run_ab_stripping.m：全流程调度。
mBDS_B1I、mBDS-3_B1C：两套软件接收机及码生成、跟踪、译码依赖。
convert_labsat3w.m：LS3W转换；tools/receiverSegment.m：强信号接收机调用。
tools/alignAidsToB.m：对齐与辅助符号转移。
LS3W_Real_B1I_B1C/Joint_B1I_B1C_PhaseComp_Acquisition.m：联合捕获入口。
tools/aidedJointAcquisition.m：辅助剥离、相位补偿、加权检测。
其余tools文件是导出、校验、统计、加噪及绘图依赖。
run_labsat_messages.m仅处理强信号电文；validate_cn0_attenuation.m用于
加噪验证；compare_b1i_b1c.m用于两路载噪比分析，均为可选入口。

6. 输出位置
output/AB_File048_demo/A：A电文和覆盖码TXT。
output/AB_File048_demo/B_prefix：用于对齐的B强前缀符号TXT。
output/AB_File048_demo/B_aligned：覆盖B、供剥离使用的辅助TXT。
reports：A载噪比等统计；results/AB_acquisition：捕获CSV、PNG、FIG。
cache：自动生成的IQ、跟踪、对齐和搜索缓存；键一致时可以复用。
output只放电文及覆盖码。再次公开时不要附带运行生成的数据目录。

同一配置完成全流程后，只重跑联合捕获：
  addpath(fullfile(pwd,'LS3W_Real_B1I_B1C'));
  r = Joint_B1I_B1C_PhaseComp_Acquisition();
该调用重新读取配置文件；前一轮临时改过cfg时，须先把同样的B路径、
起点、长度和PRN等保存到配置文件。移动录制后应重新生成绑定来源的辅助TXT。
