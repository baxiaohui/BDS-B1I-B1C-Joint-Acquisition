function [eph, subFrameStart, SOW, messages] = NAVdecoding(prompt,PRN)
% All complete navigation pages; retain the original three outputs.
% CU Multi-GNSS SDR; original NAV decoding by Yafeng Li, Nagaraj C.
% Shivaramaiah and Dennis M. Akos. GPL version 2 or later; no warranty.
eph=eph_structure_init(); SOW=inf;
[messages,subFrameStart]=decodeB1IMessages(prompt,PRN);
if isempty(messages), return; end
SOW=messages(1).sow;
% Legacy ephemeris parser requires 30 s: 5 D1 subframes / 50 D2 pages.
if PRN<=5 || PRN>=59, count=50; else, count=5; end
if numel(messages)<count || ~all([messages(1:count).valid]), return; end
[eph,SOW]=ephemeris([messages(1:count).rawBits],PRN);
SOW=messages(1).sow;
end
