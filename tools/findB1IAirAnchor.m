function anchor = findB1IAirAnchor(t,s)
% Read only a preamble + first two words (60 bits), not a complete frame.
o=demodB1IAir(t,s); bits=o.data(:,2).'; pre=[1 1 1 0 0 0 1 0 0 1 0];
anchor=[];
for invert=0:1
    b=xor(bits,invert);
    corr=conv(2*double(b)-1,fliplr(2*pre-1),'valid');
    for k=find(corr==11)
        if k+59>numel(b), continue; end
        w=double(b(k+(0:59))); codes=zeros(3,15); codes(1,:)=w(16:30);
        word=w(31:60);
        codes(2,:)=[word(1:2:22) word(23:2:30)];
        codes(3,:)=[word(2:2:22) word(24:2:30)];
        [dec,err]=bchdec(gf(codes,1),15,11);
        if any(err~=0), continue; end % strict short-header anchor
        payload=[w(1:15) reshape(double(dec.x).',1,[])];
        id=bin2dec(char(payload(16:18)+48));
        sow=bin2dec(char(payload(19:38)+48));
        if id<1 || id>5 || sow>=604800, continue; end
        if t.PRN<=5 || t.PRN>=59, sow=mod(sow+0.6*(id-1),604800); end
        anchor=struct('startIndex',o.startIndices(k),'sow',sow,'week',nan, ...
            'rawBits',zeros(1,300));
        return;
    end
end
end
