close all

Nbits = 10;

bits1 = randi([0 1], 1, Nbits);
bits2 = randi([0 1], 1, Nbits);

sign1 = bits1 * 2 - 1;
sign2 = bits2 * 2 - 1;

code1 = [1 1 1 1 1 1 1 1];
code2 = [1 -1 1 -1 1 -1 1 -1];
%code2 = [1 -1 -1 1 -1 1 1 -1];

l = length(code1);

out1 = [];
out2 = [];

for i=1:Nbits
    out1 = [out1 code1*sign1(i)];
    out2 = [out2 code2*sign2(i)];
end

out = out1+out2;

SNR = 30;

in = awgn(out, SNR);

in = [zeros(1,l) in zeros(1,l)];

r1 = [];
r2 = [];

for j=1:(Nbits+1)*l
    
    sum1 = 0;
    sum2 = 0;
    
    for i=1:l
        sum1 = sum1 + in(j+i)*code1(i);
        sum2 = sum2 + in(j+i)*code2(i);
    end
    
    r1 = [r1 sum1];
    r2 = [r2 sum2];
end

x = zeros(1, (Nbits+1)*l);
for i=1:(Nbits+1)*l
    x(i)=i;
end

figure
stem(x,r1)
figure
stem(x,r2)

data1 = [];
data2 = [];

for i=1:length(r1)
    if r1(i) > 5
        data1 = [data1 1];
    elseif r1(i) < -5
        data1 = [data1 0];
    end
end

for i=1:length(r2)
    if r2(i) > 5
        data2 = [data2 1];
    elseif r2(i) < -5
        data2 = [data2 0];
    end
end

er1 = 0;
er_xor1 = xor(bits1, data1);

er2 = 0;
er_xor2 = xor(bits2, data2);

for i=1:length(bits1)
    er1 = er1 + er_xor1(i);
    er2 = er2 + er_xor2(i);
end

er1
er2