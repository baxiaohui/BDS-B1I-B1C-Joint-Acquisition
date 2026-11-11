function observed = demodB1IAir(t,s)
% NH/bit synchronization and raw on-air decisions; NO frame/FEC required.
z=complex(double(t.I_P(:)),double(t.Q_P(:))); n=numel(z);
rx=t.absoluteSample(:)/s.samplingFreq-t.remCodePhase(:)./t.codeFreq(:) ...
    -s.filterDelaySeconds-s.segmentStartSecond;
geo=t.PRN<=5 || t.PRN>=59;
if geo, code=1; stride=2; phases=0:1;
else, code=1-2*('00000100110101001110'-48); stride=20; phases=0:19; end
scores=zeros(size(phases));
warmup=min(3000,floor(n/4));
for h=phases
    first=mod(-h,stride)+1; starts=first:stride:n-stride+1;
    starts=starts(starts>warmup); ix=(0:stride-1)'+starts;
    blocks=z(ix);
    if ~geo, blocks=blocks.*code(:); end
    scores(h+1)=mean(abs(sum(blocks,1)).^2);
end
[peak,which]=max(scores); h0=phases(which);
if ~geo
    assert(peak/max(median(scores),eps)>2,'C%02d: NH synchronization is ambiguous.',t.PRN);
end
first=mod(-h0,stride)+1; starts=first:stride:n-stride+1;
ix=(0:stride-1)'+starts; blocks=z(ix);
if ~geo, blocks=blocks.*code(:); end
soft=sum(blocks,1); rotation=exp(-0.5j*angle(sum(soft.^2)));
bits=double(real(soft*rotation)<0);
ends=interp1(1:n,rx,starts+stride,'linear','extrap');
observed.data=[nan(numel(starts),1) bits(:) 1-2*bits(:) rx(starts) ...
    ones(numel(starts),1) ends(:) zeros(numel(starts),1)];
phase=mod(h0+(0:n-1),numel(code))+1; pm=code(phase);
overlayEnd=interp1(1:n,rx,2:n+1,'linear','extrap');
observed.overlay=[nan(n,1) (1-pm(:))/2 pm(:) rx phase(:) overlayEnd(:) ones(n,1)];
observed.startIndices=starts; observed.step=stride*s.intTime;
observed.nhIndexAtFirstEpoch=h0+1;
end
