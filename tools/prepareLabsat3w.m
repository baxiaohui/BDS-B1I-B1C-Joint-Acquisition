function prepared = prepareLabsat3w(fileName, cacheDir, seconds, iqSign)
% Extract RF channel A, downconvert and FIR-decimate into interleaved int8.
% B1I: 5.8 Msps; B1C: 14.5 Msps, BOC(1,1) narrowband tracking.
% iqSign=-1 conjugates the input to test a reversed recording convention.
if nargin < 4, iqSign = 1; end
assert(ismember(iqSign,[-1 1]), 'iqSign must be +1 or -1.');
info = readLabsat3wInfo(fileName);
assert(info.SMP == 58e6, 'This conversion profile requires 58 Msps.');
assert(seconds > 0 && seconds <= info.durationSeconds, 'Duration exceeds input.');
% Distinct recordings must never overwrite the same B1I/B1C_121s files.
% Include source identity/format and conversion version in the namespace.
digest=java.security.MessageDigest.getInstance('SHA-256');
identity=struct('schema',2,'source',info,'profile','A_fir320_kaiser7_cut2.3_5.8');
digest.update(unicode2native(jsonencode(identity),'UTF-8'));
id=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
[~,stem]=fileparts(fileName);
cacheDir=fullfile(cacheDir,[stem '_' id(1:12)]);
if ~isfolder(cacheDir), mkdir(cacheDir); end
prepared = struct('signal',{},'fileName',{},'samplingFreq',{},'seconds',{}, ...
    'source',{},'iqSign',{},'filterDelaySeconds',{});
signals = {'B1I','B1C'};
fRF = [1561.098e6 1575.42e6];
factors = [10 4];
cutoffs = [2.3e6 5.8e6];
decimators = cell(1,2);
outputs = zeros(1,2);
for k = 1:2
    name = fullfile(cacheDir,sprintf('%s_%gs_iq%+d.bin',signals{k},seconds,iqSign));
    prepared(k) = struct('signal',signals{k},'fileName',name, ...
        'samplingFreq',info.SMP/factors(k),'seconds',seconds,'source',info, ...
        'iqSign',iqSign,'filterDelaySeconds',160/info.SMP);
end
meta = fullfile(cacheDir,sprintf('prepared_%gs_iq%+d.mat',seconds,iqSign));
if isfile(meta)
    old = load(meta,'prepared','cacheSignature');
    good = isequal(old.prepared,prepared) && isfield(old,'cacheSignature');
    for k=1:2
        d=dir(prepared(k).fileName);
        good=good && ~isempty(d) && d.bytes==2*floor(seconds*prepared(k).samplingFreq);
        if good
            good=d.bytes==old.cacheSignature(k).bytes && ...
                d.datenum==old.cacheSignature(k).modifiedDatenum;
        end
    end
    if good, fprintf('Using verified source-isolated conversion cache.\n'); return; end
end
% Invalidate a stale completion marker before any replacement. Write to
% private staging files; interrupted conversions cannot look complete.
if isfile(meta), delete(meta); end
staged={ [tempname(cacheDir) '.bin.part'],[tempname(cacheDir) '.bin.part'],[tempname(cacheDir) '.mat'] };
stageCleanup=onCleanup(@()removeStaging(staged));
for k = 1:2
    b = fir1(320,cutoffs(k)/(info.SMP/2),kaiser(321,7));
    decimators{k} = dsp.FIRDecimator(factors(k),b);
    [outputs(k),msg] = fopen(staged{k},'wb');
    if outputs(k)<=0, closeOutputs(outputs); error('Cannot open staged IQ: %s',msg); end
end
outCleanup=onCleanup(@() closeOutputs(outputs));
[fid,msg] = fopen(fileName,'rb','ieee-le');
assert(fid>0,'%s',msg);
inCleanup=onCleanup(@() fclose(fid));
totalRegisters = floor(seconds*info.SMP/10);
chunkRegisters = 58000; % 10 ms, bounded memory, divisible by both factors
done = 0;
while done < totalRegisters
    count = min(chunkRegisters,totalRegisters-done);
    regs = fread(fid,count,'*uint64');
    assert(numel(regs)==count,'Unexpected LS3W EOF.');
    x = unpackLabsat3wA(regs);
    if iqSign < 0, x=conj(x); end
    n = double(done*10)+(0:numel(x)-1).';
    for k = 1:2
        freq = fRF(k)-info.CFA;
        carrier = exp(complex(zeros(size(n),'single'), ...
            single(-2*pi*mod(n*freq/info.SMP,1))));
        y = decimators{k}(x.*carrier);
        iq = zeros(2,numel(y),'int8');
        iq(1,:)=int8(round(32*real(y)));
        iq(2,:)=int8(round(32*imag(y)));
        assert(fwrite(outputs(k),iq,'int8')==numel(iq),'Failed to write IQ cache.');
    end
    done=done+count;
    if mod(done,chunkRegisters*500)==0 || done==totalRegisters
        fprintf('LS3W conversion: %.2f / %.2f s\n',done*10/info.SMP,seconds);
    end
end
clear outCleanup inCleanup;
assert(isequal(info,readLabsat3wInfo(fileName)),'LS3W source changed during conversion.');
cacheSignature=repmat(struct('bytes',0,'modifiedDatenum',0),1,2);
for k=1:2
    [ok,msg]=movefile(staged{k},prepared(k).fileName,'f'); assert(ok,'%s',msg);
    d=dir(prepared(k).fileName);
    assert(d.bytes==2*floor(seconds*prepared(k).samplingFreq),'Incomplete staged IQ conversion.');
    cacheSignature(k)=struct('bytes',d.bytes,'modifiedDatenum',d.datenum);
end
save(staged{3},'prepared','cacheSignature');
[ok,msg]=movefile(staged{3},meta,'f'); assert(ok,'%s',msg);
clear stageCleanup;
end

function closeOutputs(fids)
for fid=fids
    if fid>0 && ~isempty(fopen(fid)), fclose(fid); end
end
end
function removeStaging(files)
% Only these generated staging files, never final IQ or a source recording.
for k=1:numel(files)
    if isfile(files{k}), delete(files{k}); end
end
end
