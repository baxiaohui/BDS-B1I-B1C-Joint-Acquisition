function [ok, remainder] = crc24qCheck(bits)
% CRC-24Q, MSB first, initial state zero, no final XOR; includes CRC bits.
state=uint32(0); polynomial=uint32(hex2dec('864CFB'));
for bit=double(bits(:).')
    feedback=xor(bitget(state,24)~=0,bit~=0);
    state=bitand(bitshift(state,1),uint32(hex2dec('FFFFFF')));
    if feedback, state=bitxor(state,polynomial); end
end
remainder=state; ok=state==0;
end
