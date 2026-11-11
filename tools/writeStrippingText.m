function files = writeStrippingText(track,messages,settings,signal,outputDir)
% Simple parallel timestamped files for navigation/overlay wipe-off.
if isfield(settings,'exportAirBits') && settings.exportAirBits
    files=writeAirStrippingText(track,messages,settings,signal,outputDir);
    return;
end
valid=messages([messages.valid]);
if isempty(valid)
    files=struct('dataBits','','overlay','');
    return;
end
base=sprintf('%s_PRN%02d',signal,track.PRN);
extended=isfield(settings,'exportReceiveTime') && settings.exportReceiveTime;
if extended
    % Remove the fractional-chip remainder as well as the FIR group delay.
    rx=double(track.absoluteSample)/settings.samplingFreq - ...
        double(track.remCodePhase)./double(track.codeFreq)-settings.filterDelaySeconds;
    if isfield(settings,'segmentStartSecond'), rx=rx-settings.segmentStartSecond; end
end
files.dataBits=fullfile(outputDir,[base '_data_bits.txt']);
if strcmp(signal,'B1I')
    files.overlay=fullfile(outputDir,[base '_NH.txt']);
    geo=track.PRN<=5 || track.PRN>=59;
    if geo, code=1; dt=0.002; kind='D2_NO_NH (constant +1)';
    else
        code=1-2*('00000100110101001110'-48);
        dt=0.020; kind='D1_NH (1 ms/chip, 20 chips/period)';
    end
else
    files.overlay=fullfile(outputDir,[base '_secondary_code.txt']);
    code=generate2ndCode(track.PRN); dt=0.010;
    kind='B1C_PILOT_SECONDARY (10 ms/chip, 1800 chips/period)';
end
[fid,msg]=fopen(files.dataBits,'wt','n','UTF-8'); assert(fid>0,'%s',msg);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# %s PRN%02d DATA; re-encoded from validated messages; NO NH/secondary multiplication\n',signal,track.PRN);
fprintf(fid,'# bit_duration_seconds=%.3f; PM1=1-2*BIT (0 -> +1, 1 -> -1)\n',dt);
fprintf(fid,'# BDT_SOW is transmit time at the start of this on-air bit/symbol.\n');
fprintf(fid,'# valid_frames=%d; gaps are not filled with predicted data bits.\n',numel(valid));
writeWeek(fid,valid);
if extended
    fprintf(fid,'# RX_SECONDS is relative to this processing window; may start before zero for a spanning bit.\n');
    fprintf(fid,'BDT_SOW\tBIT\tPM1\tRX_SECONDS\tPHASE_INDEX\n');
else, fprintf(fid,'BDT_SOW\tBIT\tPM1\n'); end
for m=valid
    raw=reencodeNavMessage(signal,m.decodedBits);
    times=mod(m.sow+(0:numel(raw)-1)*dt,604800);
    if extended
        stride=round(dt/settings.intTime);
        idx=m.startIndex+(0:numel(raw)-1)*stride;
        fprintf(fid,'%.3f\t%d\t%d\t%.9f\t%d\n', ...
            [times;raw;1-2*raw;rx(idx);1:numel(raw)]);
    else, fprintf(fid,'%.3f\t%d\t%d\n',[times;raw;1-2*raw]); end
end
clear cleanup;
[fid,msg]=fopen(files.overlay,'wt','n','UTF-8'); assert(fid>0,'%s',msg);
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# %s PRN%02d %s\n',signal,track.PRN,kind);
fprintf(fid,'# PM1=1-2*BIT; multiply matching-time correlator samples by PM1 to wipe off this code.\n');
fprintf(fid,'# BDT_SOW is transmit time; phase is derived from a validated frame anchor.\n');
writeWeek(fid,valid);
anchor=valid(1); step=settings.intTime;
% Confirm that all validated anchors agree with the continuous epoch count.
for m=valid
    expected=mod(anchor.sow+(m.startIndex-anchor.startIndex)*step,604800);
    delta=mod(m.sow-expected+302400,604800)-302400;
    assert(abs(delta)<step/2,'Inconsistent frame time / code phase; do not export overlay.');
end
fprintf(fid,'# Anchor BDT_SOW=%.3f at tracking index %d.\n',anchor.sow,anchor.startIndex);
fprintf(fid,'# Earlier/later chips are extrapolated over continuous tracking, not independent lock checks.\n');
if extended, fprintf(fid,'BDT_SOW\tBIT\tPM1\tRX_SECONDS\tPHASE_INDEX\n');
else, fprintf(fid,'BDT_SOW\tBIT\tPM1\n'); end
n=numel(track.I_P); indices=1:n;
phase=mod(indices-anchor.startIndex,numel(code))+1;
pm1=double(code(phase));
times=mod(anchor.sow+(indices-anchor.startIndex)*step,604800);
if extended
    fprintf(fid,'%.3f\t%d\t%d\t%.9f\t%d\n',[times;(1-pm1)/2;pm1;rx;phase]);
else, fprintf(fid,'%.3f\t%d\t%d\n',[times;(1-pm1)/2;pm1]); end
end

function writeWeek(fid,valid)
if ~isempty(valid) && isfield(valid,'week')
    weeks=[valid.week]; weeks=weeks(isfinite(weeks));
    if ~isempty(weeks), fprintf(fid,'# BDT_WEEK=%d\n',weeks(1)); end
end
end
