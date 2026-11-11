function settings = prepareReceiverInput(settings)
% Adapt an LS3W setting before any legacy binary fread/fseek is called.
% Calling again after conversion is safe: fileName already points to IQ BIN.
[~,~,extension]=fileparts(settings.fileName);
if ~strcmpi(extension,'.LS3W'), return; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root);
assert(isfield(settings,'signal') && ismember(string(settings.signal),["B1I","B1C"]), ...
    'Set settings.signal to B1I or B1C for LS3W input.');
assert(settings.msToProcess>0 && mod(settings.msToProcess,10)==0, ...
    'msToProcess must be positive and divisible by 10.');
assert(settings.skipNumberOfSamples>=0 && settings.skipNumberOfSamples==fix(settings.skipNumberOfSamples), ...
    'skipNumberOfSamples must be a nonnegative integer of converted IQ samples.');
if ~isfield(settings,'labsatCacheDir'), settings.labsatCacheDir=fullfile(root,'cache'); end
if strcmp(settings.signal,'B1I'), outFs=5.8e6; else, outFs=14.5e6; end
% Leave one full second for the legacy parallel tracker's shared read buffer.
seconds=ceil((settings.msToProcess/1000+settings.skipNumberOfSamples/outFs+1)*1000)/1000;
converted=convert_labsat3w(settings.fileName,seconds,settings.labsatCacheDir);
input=converted.(settings.signal);
settings.originalFileName=input.source.fileName;
settings.originalSamplingFreq=input.source.SMP;
settings.sourceInfo=input.source;
settings.filterDelaySeconds=input.filterDelaySeconds;
settings.fileName=input.fileName;
settings.dataType=input.dataType;
settings.fileType=input.fileType;
settings.samplingFreq=input.samplingFreq;
settings.IF=input.IF;
settings.inputFormat='IQ_BIN';
if strcmp(settings.signal,'B1C')
    % The decimation filter retains BOC(1,1); do not enable BOC(6,1) tracking.
    settings.fullBandEn=0; settings.FEBW=11.6e6;
end
fprintf('%s input adapted from LS3W to receiver-ready IQ BIN.\n',settings.signal);
end
