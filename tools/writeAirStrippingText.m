function files = writeAirStrippingText(track,messages,settings,signal,out)
% Output EVERY observed on-air bit, independent of complete-frame coverage.
% Validated/re-encoded symbols replace hard decisions where available.
valid=messages([messages.valid]);
files=struct('dataBits','','overlay','');
if ~isfield(settings,'filterDelaySeconds'), settings.filterDelaySeconds=0; end
if ~isfield(settings,'segmentStartSecond'), settings.segmentStartSecond=0; end
n=numel(track.I_P); step=settings.intTime;
rx=track.absoluteSample/settings.samplingFreq-track.remCodePhase./track.codeFreq ...
    -settings.filterDelaySeconds-settings.segmentStartSecond;
if strcmp(signal,'B1I')
    if track.PRN<=5 || track.PRN>=59, code=1; stride=2;
    else, code=1-2*('00000100110101001110'-48); stride=20; end
    suffix='_NH.txt';
else, code=double(generate2ndCode(track.PRN)); stride=1; suffix='_secondary_code.txt'; end
code=code(:).';
if ~isempty(valid), anchor=valid(1);
elseif strcmp(signal,'B1I')
    anchor=findB1IAirAnchor(track,settings);
else
    biFile=fullfile(fileparts(out),'B1I',sprintf('B1I_PRN%02d_NH.txt',track.PRN));
    anchor=[];
    if isfile(biFile)
        bi=readStrippingAid(biFile);
        sow=interp1(bi.rows(:,4),unwrap(bi.rows(:,1)*2*pi/604800)*604800/(2*pi),rx(1),'linear','extrap');
        sow=round(sow/0.01)*0.01;
        h=mod(round(sow/0.01),1800);
        anchor=struct('startIndex',1-h,'sow',mod(sow-h*0.01,604800), ...
            'week',bi.week,'rawBits',zeros(1,1800));
        phase=mod((1:n)-anchor.startIndex,1800)+1;
        [~,phaseCheck]=resolveB1CCarrierPolarity( ...
            complex(track.I_P,track.Q_P),complex(track.Pilot_I_P,track.Pilot_Q_P), ...
            code(phase),rx,step);
        if ~phaseCheck.valid, anchor=[]; end
    end
end
if isempty(anchor), return; end % No time/code lock: never fabricate phases.
code=code(:).';
phase=mod((1:n)-anchor.startIndex,numel(code))+1;
pm=code(phase); tx=mod(anchor.sow+((1:n)-anchor.startIndex)*step,604800);
base=sprintf('%s_PRN%02d',signal,track.PRN);
files.dataBits=fullfile(out,[base '_data_bits.txt']); files.overlay=fullfile(out,[base suffix]);
weeks=[valid.week]; weeks=weeks(isfinite(weeks)); week=anchor.week;
if ~isempty(weeks), week=weeks(1); end
% Include the partially observed first/last bit. A partial bit's decision
% uses only available prompts; it is explicitly marked FEC_OK=0.
first=anchor.startIndex+floor((1-anchor.startIndex)/stride)*stride;
starts=first:stride:n; bits=zeros(size(starts)); fec=zeros(size(starts));
soft=complex(track.I_P,track.Q_P);
if strcmp(signal,'B1I'), soft=soft.*pm; end
summed=complex(zeros(size(starts)));
for j=1:numel(starts)
    ix=max(1,starts(j)):min(n,starts(j)+stride-1); summed(j)=sum(soft(ix));
end
if strcmp(signal,'B1C')
    % Pilot-assisted phase reference also resolves carrier pi flips inside
    % a recording. A single rotation over the whole recording cannot do so.
    [recovered,phaseCheck]=resolveB1CCarrierPolarity( ...
        complex(track.I_P,track.Q_P),complex(track.Pilot_I_P,track.Pilot_Q_P),pm,rx,step);
    assert(phaseCheck.valid, ...
        'B1C PRN%02d secondary/phase verification failed before air-bit export (block %.3f, phase %.3f).', ...
        track.PRN,phaseCheck.coherence,phaseCheck.dataPilotSquaredCoherence);
    bits=recovered.';
else
    % B1I keeps the existing NH/bit integration and global polarity step.
    rotation=exp(-0.5j*angle(sum(summed.^2))); bits=double(real(summed*rotation)<0);
end
% Resolve the global bit polarity with already validated on-air symbols.
expected=[]; observed=[];
for m=valid
    raw=reencodeNavMessage(signal,m.decodedBits);
    j=round((m.startIndex-first)/stride)+1+(0:numel(raw)-1);
    expected=[expected raw]; observed=[observed bits(j)]; %#ok<AGROW>
end
if ~isempty(expected) && sum(expected~=observed)>numel(expected)/2, bits=1-bits; end
for m=valid
    raw=reencodeNavMessage(signal,m.decodedBits);
    j=round((m.startIndex-first)/stride)+1+(0:numel(raw)-1);
    bits(j)=raw; fec(j)=1;
end
bitTx=mod(anchor.sow+(starts-anchor.startIndex)*step,604800);
bitRx=interp1(1:n,rx,starts,'linear','extrap');
bitEnd=interp1(1:n,rx,starts+stride,'linear','extrap');
rows=[bitTx;bits;1-2*bits;bitRx;mod(round((starts-anchor.startIndex)/stride), ...
    numel(anchor.rawBits))+1;bitEnd;fec].';
writeRows(files.dataBits,rows,signal,track.PRN,week,stride*step,'AIR_DATA',anchor);
ends=interp1(1:n,rx,2:n+1,'linear','extrap');
overlay=[tx;(1-pm)/2;pm;rx;phase;ends;ones(1,n)].';
writeRows(files.overlay,overlay,signal,track.PRN,week,step,'OVERLAY',anchor);
end
function writeRows(file,rows,signal,prn,week,step,kind,anchor)
fid=fopen(file,'wt','n','UTF-8'); assert(fid>0,'Cannot write %s.',file);
closeFile=onCleanup(@()fclose(fid));
fprintf(fid,'# %s PRN%02d %s; all observed on-air bits/symbols\n',signal,prn,kind);
fprintf(fid,'# BDT_WEEK=%g; bit_duration_seconds=%.3f; PM1=1-2*BIT\n',week,step);
fprintf(fid,'# Anchor BDT_SOW=%.3f at tracking index %d.\n',anchor.sow,anchor.startIndex);
fprintf(fid,'# FEC_OK=1: checked/re-encoded frame; 0: hard decision, including partial edge frames.\n');
fprintf(fid,'# For OVERLAY, FEC_OK=1 denotes synchronized deterministic code phase.\n');
fprintf(fid,'# Frame validation is NOT an output gate. RX_SECONDS is relative to this processing window.\n');
fprintf(fid,'BDT_SOW\tBIT\tPM1\tRX_SECONDS\tPHASE_INDEX\tEND_SECONDS\tFEC_OK\n');
fprintf(fid,'%.3f\t%d\t%d\t%.9f\t%d\t%.9f\t%d\n',rows.');
end
