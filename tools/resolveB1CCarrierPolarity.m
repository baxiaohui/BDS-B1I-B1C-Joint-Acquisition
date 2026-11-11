function [bits,check] = resolveB1CCarrierPolarity(data,pilot,secondary,rx,intTime,use)
%RESOLVEB1CCARRIERPOLARITY Verify secondary code and resolve carrier pi flips.
% Verify known secondary signs over one-second blocks BEFORE using the
% wiped pilot as a local phase reference. This never fits navigation bits.
data=data(:); pilot=pilot(:); secondary=secondary(:); rx=rx(:);
n=numel(data);
assert(numel(pilot)==n && numel(secondary)==n && numel(rx)==n && n>=10, ...
    'B1C data/pilot/secondary/time lengths must agree.');
assert(all(isfinite(data)) && all(isfinite(pilot)) && all(isfinite(rx)) && ...
    all(diff(rx)>0) && all(ismember(secondary,[-1 1])) && intTime>0, ...
    'Invalid B1C phase-verification input.');
if nargin<6
    warmup=min(3,max(0,(rx(end)-rx(1))/4)); use=rx>=rx(1)+warmup;
end
use=logical(use(:)); assert(numel(use)==n && nnz(use)>=20,'Too few B1C phase-check symbols.');
p=pilot.*secondary;
groups=floor(rx(use)); uniqueGroups=unique(groups);
blockSum=complex(zeros(numel(uniqueGroups),1)); blockPower=zeros(size(blockSum));
counts=zeros(numel(uniqueGroups),1); selected=find(use);
for k=1:numel(uniqueGroups)
    ix=selected(groups==uniqueGroups(k));
    blockSum(k)=sum(p(ix)); blockPower(k)=sum(abs(p(ix))); counts(k)=numel(ix);
end
complete=counts>=max(20,round(.5/intTime));
score=abs(blockSum)./max(real(blockPower),eps);
coherence=sum(abs(blockSum(complete)))/max(sum(real(blockPower(complete))),eps);
coverage=sum(counts(complete))/nnz(use);
goodFraction=sum(counts(complete & score>=.75))/max(sum(counts(complete)),1);
% Squaring removes BPSK data signs and carrier pi ambiguities. This finds
% data-minus-pilot phase modulo pi; one GLOBAL polarity is resolved later
% against independently validated A air bits / FEC frames.
product=data(use).^2.*conj(p(use).^2);
coupling=abs(sum(product))/max(sum(abs(product)),eps);
relativePhase=.5*angle(sum(product));
window=2*floor((.05/intTime)/2)+1; % odd, about 50 ms; same raw samples
phasePilot=movsum(p,window);
rotation=exp(-1j*(angle(phasePilot)+relativePhase));
bits=double(real(data.*rotation)<0);
jumps=abs(angle(phasePilot(2:end).*conj(phasePilot(1:end-1))))>pi/2;
check=struct('valid',coherence>=.75 && goodFraction>=.8 && coverage>=.8 && coupling>=.5, ...
    'coherence',coherence,'globalCoherence',abs(sum(p(use)))/max(sum(abs(p(use))),eps), ...
    'goodBlockFraction',goodFraction,'checkedSymbolFraction',coverage, ...
    'dataPilotSquaredCoherence',coupling,'relativeDataPilotPhaseRad',relativePhase, ...
    'phaseSmoothingSymbols',window,'phasePolarityTransitions',sum(jumps), ...
    'blockStartSecond',uniqueGroups,'blockCoherence',score,'blockSymbolCount',counts);
end
