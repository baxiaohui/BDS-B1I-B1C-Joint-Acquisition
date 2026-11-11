function reports = run_labsat_messages(outputDir, options)
%RUN_LABSAT_MESSAGES Extract B1I/B1C navigation messages from File_048.LS3W.
%   run_labsat_messages
%   run_labsat_messages('D:\my_navigation_txt')
%   options fields: msToProcess (45000), maxChannels (4), signal ('both'),
%   prns ([] means all), inputFile (File_048.LS3W), reuseTracking (true).
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(outputDir), outputDir=fullfile(root,'output'); end
if nargin<2, options=struct(); end
defaults=struct('msToProcess',45000,'maxChannels',4,'signal','both', ...
    'prns',[],'inputFile',fullfile(root,'data','File_048.LS3W'),'reuseTracking',true);
for f=fieldnames(defaults).'
    if ~isfield(options,f{1}), options.(f{1})=defaults.(f{1}); end
end
assert(options.msToProcess>=1000 && mod(options.msToProcess,10)==0, ...
    'msToProcess must be at least 1000 and divisible by 10.');
assert(ismember(lower(string(options.signal)),["both","b1i","b1c"]),'Invalid signal.');
assert(options.maxChannels>=1 && options.maxChannels==fix(options.maxChannels),'Invalid channel count.');
oldPath=path; cleanup=onCleanup(@()path(oldPath));
addpath(fullfile(root,'tools'));
runId=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
% One extra second covers the shared tracking read buffer and code offset.
seconds=options.msToProcess/1000+1;
converted=convert_labsat3w(options.inputFile,seconds,fullfile(root,'cache'));
prepared=[converted.B1I converted.B1C];
reports=struct();
for s=1:2
    signal=prepared(s).signal;
    if ~strcmpi(options.signal,'both') && ~strcmpi(options.signal,signal), continue; end
    if s==1, receiver='mBDS_B1I'; else, receiver='mBDS-3_B1C'; end
    path(oldPath); addpath(fullfile(root,'tools'));
    addpath(fullfile(root,receiver),fullfile(root,receiver,'include'), ...
        fullfile(root,receiver,'Common'));
    clear initSettings preRun trkChannelsParallel eph_structure_init ephemeris;
    clear corrMatlabParallelBPSK corrMatlabParallelB1C;
    settings=initSettings();
    settings.fileName=prepared(s).fileName;
    settings.inputFormat='IQ_BIN';
    settings.dataType='int8'; settings.fileType=2;
    settings.samplingFreq=prepared(s).samplingFreq; settings.IF=0;
    settings.correlatorType=0; settings.trkMode=1;
    settings.fullBandEn=0; settings.FEBW=11.6e6;
    if s==2, settings.dllCorrelatorSpacing=0.15; end
    settings.msToProcess=options.msToProcess;
    settings.numberOfChannels=options.maxChannels;
    settings.plotTracking=0;
    settings.showWaitbar=false;
    settings.messageOutputDir=char(outputDir); settings.messageRunId=runId;
    settings.filterDelaySeconds=prepared(s).filterDelaySeconds;
    settings.originalFileName=prepared(s).source.fileName;
    settings.originalSamplingFreq=prepared(s).source.SMP;
    settings.sourceInfo=prepared(s).source;
    acquisitionFile=fullfile(root,'cache',[signal '_acquisition.mat']);
    acqCached=false;
    if isfile(acquisitionFile)
        prior=load(acquisitionFile,'acqResults','settings');
        acqCached=prior.settings.samplingFreq==settings.samplingFreq && ...
            prior.settings.IF==settings.IF && isfield(prior.settings,'sourceInfo') && ...
            isequal(prior.settings.sourceInfo,settings.sourceInfo);
        if acqCached, acqResults=prior.acqResults; end
    end
    if ~acqCached
        fid=fopen(settings.fileName,'rb'); assert(fid>0,'Cannot open converted IQ.');
        closeInput=onCleanup(@()fclose(fid));
        if exist('gpuDeviceCount','file') && gpuDeviceCount>0, acqResults=acquisitionGPU(fid,settings);
        else, acqResults=acquisitionCPU(fid,settings); end
        clear closeInput;
        save(acquisitionFile,'acqResults','settings');
    end
    if ~isempty(options.prns)
        acqResults.carrFreq(setdiff(1:numel(acqResults.carrFreq),options.prns))=0;
    end
    channel=preRun(acqResults,settings);
    trackingFile=fullfile(root,'cache',sprintf('%s_tracking_%dms_%dch.mat', ...
        signal,options.msToProcess,options.maxChannels));
    reuse=false;
    if options.reuseTracking && isfile(trackingFile)
        prior=load(trackingFile,'trackResults','settings','channel');
        reuse=isequal(prior.channel,channel) && ...
            prior.settings.samplingFreq==settings.samplingFreq && ...
            prior.settings.msToProcess==settings.msToProcess && ...
            isequal(prior.settings.sourceInfo,settings.sourceInfo) && ...
            prior.settings.dllCorrelatorSpacing==settings.dllCorrelatorSpacing;
        if reuse, trackResults=prior.trackResults; end
    end
    if ~reuse
        showChannelStatus(channel,settings);
        fid=fopen(settings.fileName,'rb'); assert(fid>0,'Cannot open converted IQ.');
        closeInput=onCleanup(@()closeIfOpen(fid));
        tic;
        [trackResults,channel]=trkChannelsParallel(fid,channel,settings);
        clear closeInput;
        fprintf('%s tracking elapsed %.1f s\n',signal,toc);
        save(trackingFile,'trackResults','channel','settings','acqResults');
    end
    reports.(signal)=exportTrackingMessages(trackResults,settings,signal);
end
fprintf('TXT output directory: %s\n',fullfile(outputDir,runId));
end

function closeIfOpen(fid)
if ~isempty(fopen(fid)), fclose(fid); end
end
