function x = unpackLabsat3wA(registers)
% Little-endian uint64 input; four leading spare/digital bits, then ten
% chronological samples of AI,AQ,BI,BQ,CI,CQ. 0 -> +1; 1 -> -1.
% Layout checked against GNSS-SDR labsat23_source.cc (see README).
registers = registers(:).';
x = complex(zeros(10,numel(registers),'single'));
for k = 0:9
    i = 1 - 2*single(bitget(registers,60-6*k));
    q = 1 - 2*single(bitget(registers,59-6*k));
    x(k+1,:) = complex(i,q);
end
x = x(:);
end
