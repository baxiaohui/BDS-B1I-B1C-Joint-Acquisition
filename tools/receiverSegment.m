function result = receiverSegment(root, prepared, signal, segment, cfg, label)
% Run the ORIGINAL receiver acquisition/tracker/decoder on a bounded prefix.
old=path; restore=onCleanup(@()path(old));
if strcmp(signal,'B1I'), folder='mBDS_B1I'; else, folder='mBDS-3_B1C'; end
addpath(root,fullfile(root,'tools'),fullfile(root,folder), ...
    fullfile(root,folder,'include'),fullfile(root,folder,'Common'));
clear initSettings acquisitionCPU acquisitionGPU preRun trkChannelsParallel;
clear eph_structure_init ephemeris corrMatlabParallelBPSK corrMatlabParallelB1C;
settings=initSettings(); p=prepared.(signal);
settings.fileName=p.fileName; settings.inputFormat='IQ_BIN';
settings.dataType=p.dataType; settings.fileType=p.fileType;
settings.samplingFreq=p.samplingFreq; settings.IF=p.IF;
settings.originalFileName=p.source.fileName; settings.sourceInfo=p.source;
settings.filterDelaySeconds=p.filterDelaySeconds;
settings.segmentStartSecond=segment.startSecond;
settings.skipNumberOfSamples=round(segment.startSecond*p.samplingFreq);
settings.msToProcess=round(segment.processSeconds*1000);
% A needs the complete symbol spanning its nominal end to cover B up to
% that end. The caller already reserves 11 ms of real A samples. B's
% strong prefix remains an exclusive boundary: never read weak samples.
endMarginSeconds=0;
if strcmp(label,'A'), endMarginSeconds=.011; end
settings.processingEndSample=round((segment.startSecond+segment.processSeconds+endMarginSeconds)*p.samplingFreq);
settings.numberOfChannels=cfg.maxChannels;
settings.correlatorType=0; settings.trkMode=1;
settings.fullBandEn=0; settings.FEBW=11.6e6;
settings.showWaitbar=false; settings.plotTracking=0;
settings.exportReceiveTime=true;
settings.exportAirBits=true;
settings.messageOutputDir=cfg.outputDir; settings.messageRunId=label;
if ~isempty(cfg.prns), settings.acqSatelliteList=cfg.prns; end
% Broad initial search is configured separately from weak aided acquisition.
if isfield(cfg,'receiverAcquisition') && isfield(cfg.receiverAcquisition,signal)
    search=cfg.receiverAcquisition.(signal);
    assert(isscalar(search.searchBandHz) && isfinite(search.searchBandHz) && search.searchBandHz>0 && ...
        isscalar(search.frequencyStepHz) && isfinite(search.frequencyStepHz) && search.frequencyStepHz>0 && ...
        search.frequencyStepHz<=2*search.searchBandHz,'Invalid initial %s Doppler search.',signal);
    settings.acqSearchBand=search.searchBandHz;
    if strcmp(signal,'B1I')
        settings.acqSearchStep=search.frequencyStepHz;
        if isfield(search,'noncoherentMs')
            assert(isscalar(search.noncoherentMs) && isfinite(search.noncoherentMs) && ...
                search.noncoherentMs>=1 && search.noncoherentMs==fix(search.noncoherentMs), ...
                'B1I noncoherentMs must be a positive integer.');
            settings.acqNonCohTime=search.noncoherentMs;
        end
    else
        settings.acqStep=search.frequencyStepHz;
    end
end
if strcmp(signal,'B1I'), step=settings.acqSearchStep; noncoherentMs=settings.acqNonCohTime;
else, step=settings.acqStep; noncoherentMs=0; end
acquisitionOptions=struct('satellites',settings.acqSatelliteList, ...
    'searchBandHz',settings.acqSearchBand,'frequencyStepHz',step, ...
    'threshold',settings.acqThreshold,'noncoherentMs',noncoherentMs);
% Source identity is part of the cache directory; segment/PRNs part of key.
key=struct('schema',6,'signal',signal,'source',p.source, ...
    'start',segment.startSecond,'seconds',segment.processSeconds, ...
    'prns',cfg.prns,'maxChannels',cfg.maxChannels,'acquisitionOptions',acquisitionOptions, ...
    'endMarginSeconds',endMarginSeconds);
work=fullfile(cfg.cacheDir,label); if ~isfolder(work), mkdir(work); end
file=fullfile(work,[signal '_tracking.mat']);
reuse=false;
if cfg.reuseTracking && isfile(file)
    prior=load(file,'key'); reuse=isequal(prior.key,key);
end
if reuse
    data=load(file,'trackResults','channel','acqResults');
    trackResults=data.trackResults; channel=data.channel; acqResults=data.acqResults;
    fprintf('%s %s: reusing matching tracking cache.\n',label,signal);
else
    fprintf('\n%s %s: acquire/track %.3f s from file %.3f s.\n', ...
        label,signal,segment.processSeconds,segment.startSecond);
    fprintf('  PRNs=%s; Doppler=[%g,%g] Hz, step=%g Hz; threshold=%g.\n', ...
        mat2str(settings.acqSatelliteList),-settings.acqSearchBand,settings.acqSearchBand,step,settings.acqThreshold);
    fid=fopen(p.fileName,'rb'); assert(fid>0,'Cannot open IQ cache.');
    closer=onCleanup(@()closeIfOpen(fid));
    if exist('gpuDeviceCount','file') && gpuDeviceCount>0
        acqResults=acquisitionGPU(fid,settings);
    else, acqResults=acquisitionCPU(fid,settings); end
    if ~isempty(cfg.prns)
        acqResults.carrFreq(setdiff(1:numel(acqResults.carrFreq),cfg.prns))=0;
    end
    channel=preRun(acqResults,settings);
    if ~any([channel.status]=='T')
        failureFile=fullfile(work,[signal '_acquisition_failed.mat']);
        save(failureFile,'key','settings','acqResults','-v7.3');
        [peak,where]=max(acqResults.peakMetric(settings.acqSatelliteList));
        error('Receiver:NoAcquiredChannels', ...
            ['%s %s: no acquired channels in +/-%.0f Hz. Best PRN %d ratio %.4f; ' ...
             'required > %.4f. Saved diagnostic: %s'], ...
            label,signal,settings.acqSearchBand,settings.acqSatelliteList(where), ...
            peak,settings.acqThreshold,failureFile);
    end
    showChannelStatus(channel,settings);
    [trackResults,channel]=trkChannelsParallel(fid,channel,settings);
    clear closer;
    active=find([trackResults.status]=='T');
    assert(~isempty(active),'Tracking failed.');
    for j=active
        assert(all(isfinite(trackResults(j).codeFreq)) && ...
            all(diff(trackResults(j).absoluteSample)>0), ...
            'Incomplete tracking: %s PRN%d. No aid export.',signal,trackResults(j).PRN);
    end
    save(file,'key','settings','trackResults','channel','acqResults','-v7.3');
end
if strcmp(label,'B_prefix')
    % Alignment uses actual prefix prompts, with no complete-frame wait.
    report=struct('PRN',num2cell([trackResults([trackResults.status]=='T').PRN]), ...
        'valid',0,'files',struct('dataBits','','overlay',''));
else, report=exportTrackingMessages(trackResults,settings,signal); end
result=struct('report',report,'trackingFile',file,'settings',settings);
end
function closeIfOpen(fid)
if ~isempty(fopen(fid)), fclose(fid); end
end
