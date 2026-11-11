function [rxI,rxC,fs] = readLabsatWindow(file,startSecond,milliseconds,noiseModel)
% Bounded random read. Same 320th-order anti-alias filters as the receivers.
% Return both complex basebands at a COMMON 14.5 Msps; compensate FIR delay.
info=readLabsat3wInfo(file); assert(info.SMP==58e6,'Requires 58 Msps.');
fs=info.SMP/4; count=round(milliseconds*fs/1000);
first=round(startSecond*info.SMP);
assert(abs(first/info.SMP-startSecond)<1e-10,'Start must lie on a raw sample.');
rawBegin=max(0,floor((first-320)/20)*20);
rawEnd=ceil((first+(count-1)*4+161)/20)*20;
assert(rawEnd<=info.totalSamples,'Acquisition window exceeds file.');
fid=fopen(file,'rb','ieee-le'); assert(fid>0,'Cannot open LS3W.');
closer=onCleanup(@()fclose(fid));
assert(fseek(fid,rawBegin/10*8,'bof')==0,'LS3W seek failed.');
words=fread(fid,(rawEnd-rawBegin)/10,'*uint64');
assert(numel(words)==(rawEnd-rawBegin)/10,'Short LS3W read.');
x=double(unpackLabsat3wA(words));
n=rawBegin+(0:numel(x)-1).';
pick=first-rawBegin+161+(0:count-1)*4;
[filters,rf]=labsatAcquisitionFilters(info.SMP); base=cell(1,2);
for k=1:2
    carrier=exp(-1j*2*pi*mod(n*((rf(k)-info.CFA)/info.SMP),1));
    mixed=x.*carrier;
    if nargin>=4 && noiseModel.enabled
        assert(noiseModel.rawFs==info.SMP,'Noise sample-rate mismatch.');
        % Independent noise in the two disjoint RF bands, BEFORE each
        % unchanged anti-alias FIR and decimation. Do not requantize/sign().
        mixed=mixed+acquisitionNoiseSamples(rawBegin,numel(x),k,noiseModel);
    end
    y=filter(filters{k},1,mixed); base{k}=y(pick);
end
rxI=base{1}; rxC=base{2};
end
