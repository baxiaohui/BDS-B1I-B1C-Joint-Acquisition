function z = acquisitionNoiseSamples(first,count,branch,model)
% Random-access circular complex Gaussian noise. E|z|^2 = Fs * addedN0.
% Absolute raw sample indexing makes overlapping reads, FFT partition sizes,
% GPU choices and PRN ordering use EXACTLY the same noise realization.
assert(ismember(branch,[1 2]) && isscalar(branch),'branch must be 1 or 2.');
assert(first>=0 && first==fix(first) && count>=0 && count==fix(count),'Invalid sample range.');
z=complex(zeros(count,1));
if ~model.enabled || count==0 || model.addedNoisePsd(branch)==0, return; end
blockSize=65536;
stream=RandStream('Threefry','Seed',model.seed);
last=first+count-1;
for block=floor(first/blockSize):floor(last/blockSize)
    stream.Substream=2*block+branch; % 1-based substreams; branch bands independent
    a=randn(stream,2,blockSize);
    lo=max(first,block*blockSize); hi=min(last,(block+1)*blockSize-1);
    src=lo-block*blockSize+1:hi-block*blockSize+1;
    dst=lo-first+1:hi-first+1;
    z(dst)=complex(a(1,src),a(2,src)).';
end
z=z*sqrt(model.rawFs*model.addedNoisePsd(branch)/2);
% The strong prefix is never modified, even in FIR/interpolation guard reads.
z((first:first+count-1)'<model.startRawSample)=0;
end
