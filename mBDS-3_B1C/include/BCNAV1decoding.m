function [eph, firstSubFrame, TOW, messages] = BCNAV1decoding(trackResults, channelNr, settings)
% Decode every complete frame; preserve legacy ephemeris outputs.
% CU Multi-GNSS SDR; original B1C decoding by Yafeng Li, Nagaraj C.
% Shivaramaiah and Dennis M. Akos. GPL version 2 or later; no warranty.
eph=eph_structure_init(); firstSubFrame=inf; TOW=inf;
messages=struct('type',{},'startIndex',{},'rawBits',{},'decodedBits',{}, ...
    'valid',{},'status',{},'subframeId',{},'sow',{},'week',{},'ldpcOK',{});
addpath(fullfile(fileparts(mfilename('fullpath')),'ldpcDecoder'));
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'tools'));
t=trackResults(channelNr); prn=t.PRN;
secondary=double(generate2ndCode(prn)); secondary=secondary(:).';
% Narrowband pilot is in quadrature; full-band pilot is phase-rotated.
% Search both prompts so the synchronization also tolerates phase ambiguity.
pi=double(t.Pilot_I_P(:).'); pq=double(t.Pilot_Q_P(:).');
n=min([numel(pi),numel(pq),numel(t.I_P)]);
if n<1800, return; end
ci=conv(2*(pi(1:n)>0)-1,fliplr(secondary),'valid');
cq=conv(2*(pq(1:n)>0)-1,fliplr(secondary),'valid');
score=max(abs(ci),abs(cq));
candidates=find(score>=0.80*1800);
% Only one candidate per nearby correlation peak.
[~,order]=sort(score(candidates),'descend'); keep=[];
for k=candidates(order)
    if isempty(keep) || all(abs(keep-k)>10), keep(end+1)=k; end %#ok<AGROW>
end
for start=sort(keep)
    softI=t.I_P(start:start+1799); softQ=t.Q_P(start:start+1799);
    [bits,details]=B1CLDPCDecoder(softI,softQ);
    if ~isempty(bits) && bin2dec(char(bits(1:6).'+48))~=prn
        [bits,details]=B1CLDPCDecoder(-softI,-softQ);
        details.polarity=-details.polarity;
    end
    valid=false; sow=nan; week=nan; status='BCH_FAIL';
    if ~isempty(bits)
        prnOK=bin2dec(char(bits(1:6).'+48))==prn;
        sf2OK=crc24qCheck(bits(15:614)); sf3OK=crc24qCheck(bits(615:878));
        valid=prnOK && sf2OK && sf3OK;
        if valid
            week=bin2dec(char(bits(15:27).'+48));
            how=bin2dec(char(bits(28:35).'+48));
            soh=bin2dec(char(bits(7:14).'+48))*18;
            sow=how*3600+soh;
            status='CRC24Q_OK';
            eph=ephemeris(char(bits(:).'+48),eph);
            if isinf(TOW), TOW=sow; firstSubFrame=start; end
        elseif ~prnOK, status='PRN_MISMATCH';
        elseif ~sf2OK && ~sf3OK, status='CRC24Q_SF2_SF3_FAIL';
        elseif ~sf2OK, status='CRC24Q_SF2_FAIL'; else, status='CRC24Q_SF3_FAIL'; end
    end
    raw=double(softI*details.polarity>0);
    messages(end+1)=struct('type','BCNAV1','startIndex',start,'rawBits',raw, ...
        'decodedBits',double(bits(:).'),'valid',valid,'status',status, ...
        'subframeId',nan,'sow',sow,'week',week, ...
        'ldpcOK',[details.sf2SyndromeOK details.sf3SyndromeOK]); %#ok<AGROW>
end
end
