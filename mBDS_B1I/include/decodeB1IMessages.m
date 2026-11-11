function [messages, firstStart] = decodeB1IMessages(prompt, PRN)
% Synchronize every complete D1 subframe / D2 page independently of PVT.
prompt = double(prompt(:).');
assert(isscalar(PRN) && PRN>=1 && PRN<=63 && PRN==fix(PRN),'Invalid PRN.');
preamble = [1 1 1 0 0 0 1 0 0 1 0];
nh = [-1 -1 -1 -1 -1 1 -1 -1 1 1 -1 1 -1 1 -1 -1 1 1 1 -1];
geo = PRN<=5 || PRN>=59;
if geo, step=2; overlay=ones(1,2); kind='D2';
else, step=20; overlay=nh; kind='D1'; end
pattern=kron(2*preamble-1,overlay);
hard=2*(prompt>0)-1;
messages=struct('type',{},'startIndex',{},'rawBits',{},'decodedBits',{}, ...
    'valid',{},'status',{},'subframeId',{},'sow',{},'bchErrors',{},'week',{});
firstStart=inf;
if numel(prompt)<300*step+numel(pattern), return; end
corr=conv(hard,fliplr(pattern),'valid');
candidates=find(abs(corr)>=0.9*numel(pattern));
candidates=candidates(candidates>1000 & candidates+300*step<=numel(corr));
for start=candidates
    if abs(corr(start+300*step))<0.9*numel(pattern), continue; end
    b=bitDecisions(prompt(start:start+30*step-1),step,overlay);
    if isequal(b(1:11),1-preamble), b=1-b; end
    if ~isequal(b(1:11),preamble), continue; end
    [~,errors]=bchdec(gf(b(16:30),1),15,11);
    if errors>=0, firstStart=start; break; end
end
if isinf(firstStart), return; end
for start=firstStart:300*step:numel(prompt)-300*step+1
    b=bitDecisions(prompt(start:start+300*step-1),step,overlay);
    if isequal(b(1:11),1-preamble), b=1-b; end
    payload=b(1:15);
    codes=zeros(19,15); codes(1,:)=b(16:30);
    for w=1:9
        word=b(30*w+(1:30));
        codes(2*w,:)=[word(1:2:22),word(23:2:30)];
        codes(2*w+1,:)=[word(2:2:22),word(24:2:30)];
    end
    [decoded,nerr]=bchdec(gf(codes,1),15,11);
    payload=[payload reshape(double(decoded.x).',1,[])];
    id=bin2dec(char(payload(16:18)+'0'));
    sow=bin2dec(char([payload(19:26) payload(27:38)]+'0'));
    % D2 broadcasts frame-start SOW in every subframe (ICD section 5.3.3).
    if geo, sow=mod(sow+0.6*(id-1),604800); end
    week=nan;
    if id==1 && ~geo, week=bin2dec(char(payload(49:61)+'0')); end
    if id==1 && geo && bin2dec(char(payload(39:42)+'0'))==1
        week=bin2dec(char(payload(53:65)+'0'));
    end
    valid=isequal(b(1:11),preamble) && all(nerr>=0) && id>=1 && id<=5;
    if valid, status='BCH_OK'; else, status='BCH_OR_PREAMBLE_FAIL'; end
    messages(end+1)=struct('type',kind,'startIndex',start,'rawBits',double(b), ...
        'decodedBits',payload,'valid',valid,'status',status,'subframeId',id, ...
        'sow',sow,'bchErrors',double(nerr(:).'),'week',week); %#ok<AGROW>
end
end

function b=bitDecisions(p,step,overlay)
p=reshape(p,step,[]).*overlay(:);
b=double(sum(p,1)>0);
end
