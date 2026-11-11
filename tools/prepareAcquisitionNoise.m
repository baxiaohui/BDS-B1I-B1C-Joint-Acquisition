function model = prepareAcquisitionNoise(cfg)
% Effective digital C/N0 degradation, after LS3W unpacking, before FIR.
% The source recording and all tracking/alignment caches remain untouched.
cfg=resolveAcquisitionNoiseOptions(cfg); n=cfg.noise;
info=readLabsat3wInfo(cfg.ls3wFile);
[filters,~]=labsatAcquisitionFilters(info.SMP);
enbw=info.SMP*cellfun(@(b)sum(abs(b).^2),filters);
active=n.enabled && n.cn0DropDb>0;
model=struct('enabled',active,'requestedEnabled',n.enabled, ...
    'requestedDropDb',n.cn0DropDb,'appliedDropDb',double(active)*n.cn0DropDb, ...
    'seed',n.seed,'rawFs',info.SMP, ...
    'startRawSample',round((cfg.fileStartSecond+cfg.startSecond)*info.SMP), ...
    'originalNoisePsd',[NaN NaN],'addedNoisePsd',[0 0], ...
    'calibrationMode',n.calibrationMode,'bandpowerNoisePsd',[NaN NaN], ...
    'codeWeightedNoisePsd',[NaN NaN], ...
    'filterNoiseBandwidthHz',enbw,'calibrationMs',0, ...
    'calibrationPower',[NaN NaN],'calibrationMethod','disabled', ...
    'referenceCN0dBHz',n.referenceCN0dBHz,'predictedCN0dBHz',[], ...
    'domain','unpacked floating-point digital baseband, before branch FIR', ...
    'randomGenerator','Threefry; absolute raw-sample blocks; independent RF bands');
if ~isempty(n.referenceCN0dBHz)
    model.predictedCN0dBHz=n.referenceCN0dBHz-model.appliedDropDb;
end
if ~active, return; end
if isempty(n.originalNoisePsd)
    % GNSS before despreading is noise dominated. Equivalent local N0 is
    % estimated with the selected spectral weighting.
    % This includes quantization noise; it is not a measurement of RF watts.
    start=model.startRawSample/info.SMP;
    availableMs=(info.durationSeconds-start-1e-4)*1000;
    milliseconds=min(n.calibrationMs,availableMs);
    assert(milliseconds>=1,'Less than 1 ms available for noise calibration.');
    [xI,xC,fsOut]=readLabsatWindow(cfg.ls3wFile,start,milliseconds);
    power=[mean(abs(xI-mean(xI)).^2) mean(abs(xC-mean(xC)).^2)];
    model.bandpowerNoisePsd=power./enbw;
    model.codeWeightedNoisePsd=[weightedPsd(xI,fsOut,info.SMP,filters{1},1), ...
        weightedPsd(xC,fsOut,info.SMP,filters{2},2)];
    if strcmp(n.calibrationMode,'code-spectrum')
        model.originalNoisePsd=model.codeWeightedNoisePsd;
        model.calibrationMethod='Welch background PSD weighted by primary modulation spectrum and FIR response';
    else
        model.originalNoisePsd=model.bandpowerNoisePsd;
        model.calibrationMethod='noise-dominated branch power / exact FIR noise bandwidth';
    end
    model.calibrationPower=power;
    model.calibrationMs=milliseconds;
else
    model.originalNoisePsd=n.originalNoisePsd;
    model.calibrationMethod='manual originalNoisePsd';
end
assert(all(isfinite(model.originalNoisePsd)) && all(model.originalNoisePsd>0), ...
    'Cannot calibrate positive original noise PSD.');
% N0_new = N0_old * 10^(drop/10). Original noise is already in the file!
model.addedNoisePsd=model.originalNoisePsd*expm1(log(10)*n.cn0DropDb/10);
fprintf(['Acquisition-only AWGN: C/N0 drop %.6g dB, seed %u; ' ...
    'original N0 [I C]=[%.6g %.6g] sample^2/Hz (%s).\n'], ...
    n.cn0DropDb,n.seed,model.originalNoisePsd,model.calibrationMethod);
if ~isempty(model.predictedCN0dBHz)
    fprintf('  Supplied reference C/N0 [I C]=[%g %g] -> theoretical [%g %g] dB-Hz.\n', ...
        model.referenceCN0dBHz,model.predictedCN0dBHz);
end
end
function n0=weightedPsd(x,fsOut,fsRaw,b,branch)
% Noise seen by a matched correlator depends on its code spectrum, not on
% an unweighted average over the entire front-end passband. GNSS signals
% are below the pre-despreading floor; use observed PSD as background.
% Pout already contains |H|^2, hence it occurs ONLY in the denominator.
x=x-mean(x);
[p,f]=pwelch(x,hann(8192,'periodic'),4096,16384,fsOut,'centered');
H=freqz(b,1,2*pi*f/fsRaw);
if branch==1
    Tc=1/2.046e6; g=Tc*sinc(f*Tc).^2;
else
    Tc=1/1.023e6; dt=Tc/2;
    g=dt^2/Tc*sinc(f*dt).^2.*abs(1-exp(-2j*pi*f*dt)).^2;
end
n0=sum(p.*g)/sum(abs(H).^2.*g);
end
