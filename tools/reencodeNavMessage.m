function raw = reencodeNavMessage(signal, payload)
% Rebuild error-corrected on-air navigation symbols for data wipe-off.
% Input includes the CRC bits for B1C. Output excludes ranging/overlay codes.
payload=double(payload(:).');
if strcmp(signal,'B1I')
    assert(numel(payload)==224,'B1I payload must have 224 bits.');
    cw=bchenc(gf(reshape(payload(16:end),11,19).',1),15,11);
    cw=double(cw.x); raw=[payload(1:15) cw(1,:) zeros(1,270)];
    for w=1:9
        pair=cw(2*w:2*w+1,:);
        raw(w*30+(1:30))=[reshape(pair(:,1:11),1,[]) reshape(pair(:,12:15),1,[])];
    end
else
    assert(numel(payload)==878,'B1C payload must have 878 bits.');
    [~,~,h2]=getH_informationSF2(); [~,~,h3]=getH_informationSF3();
    sf2=encodeLDPC(payload(15:614),h2); sf3=encodeLDPC(payload(615:878),h3);
    data=zeros(36,48); rows3=3:3:35; rows2=setdiff(1:36,rows3);
    data(rows2,:)=reshape(sf2,48,25).'; data(rows3,:)=reshape(sf3,48,11).';
    raw=[encodeBCH(payload(1:6),21,[2 4 5 6]) ...
        encodeBCH(payload(7:14),51,[1 4 5 6 7 8]) reshape(data,1,[])];
end
end

function bits=encodeLDPC(input,H)
n=numel(input)/6;
info=bin2dec(char(reshape(input,6,[]).'+48));
rhs=gf(H(:,1:n),6)*gf(info,6);
parity=gf(H(:,n+1:end),6)\rhs;
code=[info;double(parity.x)];
assert(all(gf(H,6)*gf(code,6)==0),'LDPC re-encoding syndrome failed.');
bits=reshape((dec2bin(code,6)-48).',1,[]);
end

function bits=encodeBCH(input,count,taps)
reg=fliplr(1-2*input); bits=zeros(1,count);
for k=1:count
    bits(k)=(1-reg(end))/2;
    feedback=prod(reg(taps)); reg=circshift(reg,1); reg(1)=feedback;
end
end
