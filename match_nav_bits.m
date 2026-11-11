function matches = match_nav_bits(referenceFile, weakBits, maxErrorRate)
%MATCH_NAV_BITS Locate a short weak-signal bit sequence in a saved TXT.
% matches=match_nav_bits('.../B1I_PRN07_data_bits.txt','010101...')
% Match on-air navigation symbols using the re-encoded *_data_bits.txt.
% Positive/inverted carrier polarity are both tried.
if nargin<3, maxErrorRate=0.05; end
if ischar(weakBits) || isstring(weakBits)
    weakBits=regexprep(char(weakBits),'\s','');
    assert(all(weakBits=='0' | weakBits=='1'),'Bits must contain only 0/1.');
    weakBits=double(weakBits)-48;
else, weakBits=double(weakBits(:).'); end
assert(numel(weakBits)>=24 && all(weakBits==0 | weakBits==1), ...
    'Supply at least 24 binary bits; 60 or more is recommended.');
assert(maxErrorRate>=0 && maxErrorRate<0.5,'Invalid error rate.');
lines=readlines(referenceFile);
meta=strjoin(lines(startsWith(lines,'#')),newline);
if contains(meta,'B1C'), frameSeconds=18;
else
    prn=regexp(char(meta),'PRN[= ]?(\d+)','tokens','once');
    assert(~isempty(prn),'Cannot identify reference PRN.');
    prn=str2double(prn{1});
    if prn<=5 || prn>=59, frameSeconds=0.6; else, frameSeconds=6; end
end
lines=lines(strlength(lines)>0 & ~startsWith(lines,'#'));
assert(numel(lines)>=2,'No validated reference frames are available.');
headers=split(lines(1),sprintf('\t'));
isDetailed=any(headers=="raw_bits");
if isDetailed
    bitsColumn=find(headers=="raw_bits"); timeColumn=find(headers=="sow_bdt");
    validColumn=find(headers=="valid");
elseif any(headers=="BIT")
    bitsColumn=find(headers=="BIT"); timeColumn=find(headers=="BDT_SOW");
    token=regexp(char(meta),'bit_duration_seconds=([0-9.]+)','tokens','once');
    assert(~isempty(token),'Use the DATA bits file for navigation matching.');
    frameSeconds=str2double(token{1});
else
    bitsColumn=find(headers=="NAV_BITS"); timeColumn=find(headers=="BDT_SOW");
end
assert(~isempty(bitsColumn) && ~isempty(timeColumn),'Unsupported TXT columns.');
stream=[]; rowStart=[]; times=[];
for row=lines(2:end).'
    cols=split(row,sprintf('\t'));
    if isDetailed && cols(validColumn)~="1", continue; end
    b=char(cols(bitsColumn));
    if isempty(b), continue; end
    rowStart(end+1)=numel(stream)+1; %#ok<AGROW>
    times(end+1)=str2double(cols(timeColumn)); %#ok<AGROW>
    stream=[stream double(b)-48]; %#ok<AGROW>
end
matches=table('Size',[0 5],'VariableTypes',{'double','double','double','double','logical'}, ...
    'VariableNames',{'BDT_SOW','BitOffset','Errors','ComparedBits','Inverted'});
if numel(stream)<numel(weakBits), return; end
score=conv(2*stream-1,fliplr(2*weakBits-1),'valid');
errors=round((numel(weakBits)-abs(score))/2);
locations=find(errors<=floor(maxErrorRate*numel(weakBits)));
for loc=locations
    first=find(rowStart<=loc,1,'last');
    last=find(rowStart<=loc+numel(weakBits)-1,1,'last');
    % Do not splice across missing reference frames.
    if last>first
        gaps=mod(diff(times(first:last)),604800);
        if any(abs(gaps-frameSeconds)>1e-6), continue; end
    end
    matches(end+1,:)={times(first),loc-rowStart(first),errors(loc),numel(weakBits),score(loc)<0}; %#ok<AGROW>
end
matches=sortrows(matches,{'Errors','BDT_SOW','BitOffset'});
for k=1:height(matches)
    fprintf('BDT_SOW=%.3f  BitOffset=%d  Errors=%d/%d  Inverted=%d\n', ...
        matches.BDT_SOW(k),matches.BitOffset(k),matches.Errors(k), ...
        matches.ComparedBits(k),matches.Inverted(k));
end
fprintf('BitOffset is zero-based within the reference row; BDT_SOW identifies that row start.\n');
if height(matches)>1, fprintf('Multiple matches: use a longer sequence to disambiguate.\n'); end
end
