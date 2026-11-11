function report = exportTrackingMessages(trackResults,settings,signal)
% No minimum satellite count or positioning solution is required.
if ~isfield(settings,'messageOutputDir')
    root=fileparts(fileparts(mfilename('fullpath')));
    settings.messageOutputDir=fullfile(root,'output');
end
if ~isfield(settings,'messageRunId')
    settings.messageRunId=char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
end
out=fullfile(settings.messageOutputDir,settings.messageRunId,signal);
if ~isfolder(out), mkdir(out); end
report=struct('PRN',{},'frames',{},'valid',{},'status',{},'files',{});
for channelNr=find([trackResults.status]~='-')
    track=trackResults(channelNr);
    if strcmp(signal,'B1I'), [~,~,~,messages]=NAVdecoding(track.I_P,track.PRN);
    else, [~,~,~,messages]=BCNAV1decoding(trackResults,channelNr,settings); end
    files=writeNavText(track,messages,settings,signal,out);
    valid=sum([messages.valid]);
    if isfield(settings,'exportAirBits') && settings.exportAirBits && ~isempty(files.dataBits)
        status='AIR_BITS_SYNCED';
    elseif isempty(messages), status='NO_FRAME_SYNC';
    elseif valid==0, status='NO_VALIDATED_FRAME'; else, status='OK'; end
    if ~isempty(files.dataBits)
        fprintf('%s PRN %02d: continuous air bits (%d validated frames) -> %s\n', ...
            signal,track.PRN,valid,files.dataBits);
    else
        fprintf('%s PRN %02d: no validated frames; no output files written.\n',signal,track.PRN);
    end
    report(end+1)=struct('PRN',track.PRN,'frames',numel(messages), ...
        'valid',valid,'status',status,'files',files); %#ok<AGROW>
end
end
