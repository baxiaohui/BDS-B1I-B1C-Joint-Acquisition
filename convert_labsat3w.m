function converted = convert_labsat3w(fileName, seconds, cacheDir)
%CONVERT_LABSAT3W Convert LS3W channel A into receiver-ready B1I/B1C IQ BIN.
% converted = convert_labsat3w('File_048.LS3W',46);
% converted.B1I / converted.B1C contain fileName, dataType, fileType,
% samplingFreq, IF, source and filterDelaySeconds.
% BIN layout: signed int8, I0 Q0 I1 Q1 ... . No LS3W header or packing remains.
% Profile: 58 Msps, QUA=1, CHN=3, SFT=6; B1C uses narrowband BOC(1,1).
root=fileparts(mfilename('fullpath'));
if nargin<1 || isempty(fileName), fileName=fullfile(root,'data','File_048.LS3W'); end
if nargin<2 || isempty(seconds), seconds=46; end
if nargin<3 || isempty(cacheDir), cacheDir=fullfile(root,'cache'); end
addpath(fullfile(root,'tools'));
% Canonical paths keep the cache reusable from either receiver directory.
[ok,details]=fileattrib(fileName);
assert(ok,'LS3W input does not exist: %s',fileName);
fileName=details.Name;
prepared=prepareLabsat3w(fileName,cacheDir,seconds,1);
converted=struct();
for k=1:numel(prepared)
    item=prepared(k);
    item.dataType='int8'; item.fileType=2; item.IF=0;
    converted.(item.signal)=item;
    fprintf('%s IQ: %s | int8 I/Q | Fs=%.3f MHz | IF=0 Hz\n', ...
        item.signal,item.fileName,item.samplingFreq/1e6);
end
end
