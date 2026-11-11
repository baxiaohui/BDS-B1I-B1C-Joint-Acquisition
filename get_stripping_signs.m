function signs = get_stripping_signs(aidDir,signal,prn,bSeconds)
% Query B-aligned on-air wipe-off signs at arbitrary B-relative times.
% B1I: use dataWipePM1 (data * NH).
% B1C: use dataWipePM1 for DATA, pilotWipePM1 for PILOT, separately.
% These signs do not remove the primary ranging code or carrier.
root=fileparts(mfilename('fullpath')); addpath(fullfile(root,'tools'));
signal=upper(char(signal));
assert(ismember(string(signal),["B1I","B1C"]),'Choose B1I or B1C.');
t=bSeconds(:); assert(all(isfinite(t) & t>=0),'Invalid B-relative times.');
base=fullfile(aidDir,signal,sprintf('%s_PRN%02d',signal,prn));
data=readStrippingAid([base '_data_bits.txt']);
[signs.dataPM1,signs.dataPhaseIndex,signs.dataFecOK]=sample(data,t);
if strcmp(signal,'B1I')
    overlay=readStrippingAid([base '_NH.txt']);
    [signs.NHPM1,signs.NHPhaseIndex]=sample(overlay,t);
    signs.dataWipePM1=signs.dataPM1.*signs.NHPM1;
else
    overlay=readStrippingAid([base '_secondary_code.txt']);
    [signs.secondaryPM1,signs.secondaryPhaseIndex]=sample(overlay,t);
    signs.dataWipePM1=signs.dataPM1;
    signs.pilotWipePM1=signs.secondaryPM1;
end
end
function [pm,phase,fec]=sample(aid,t)
duration=regexp(aid.header,'B_DURATION_SECONDS=([\d.]+)','tokens','once');
if ~isempty(duration), assert(all(t<str2double(duration{1})),'Time exceeds B duration.'); end
ix=interp1(aid.rows(:,4),(1:size(aid.rows,1))',t,'previous','extrap');
assert(all(isfinite(ix)) && all(t>=aid.rows(1,4)-0.5e-9),'No aid at requested time.');
assert(all(t<aid.rows(ix,6)+0.5e-9),'Requested time falls in an unobserved aid gap.');
pm=aid.rows(ix,3); phase=aid.rows(ix,5);
if size(aid.rows,2)>=7, fec=logical(aid.rows(ix,7)); else, fec=true(size(t)); end
end
