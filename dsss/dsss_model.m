Nbits = 10;

bits = randi([0 1], 1, Nbits);

sign = bits * 2 - 1;

code = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];

out = [];

for i=1:Nbits
    out = [out code*sign(i)];
end


SNR = 30;

in = awgn(out, SNR);

in = [zeros(1,length(code)) in zeros(1,length(code))];

r = [];

for j=1:(Nbits+1)*length(code)
    
    sum = 0;
    
    for i=1:length(code)
        sum = sum + in(j+i)*code(i);
    end
    
    r = [r sum];
end

x = zeros(1, length(r));
for i=1:length(r)
    x(i)=i;
end

stem(x,r)

data = [];

for i=1:length(r)
    if r(i) > 5
        data = [data 1];
    elseif r(i) < -5
        data = [data 0];
    end
end

bits
data

er = 0;
er_xor = xor(bits, data);

for i=1:length(bits)
    er = er + er_xor(i);
end

er