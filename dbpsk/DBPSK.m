% clear all;
% close all;
% clc;

a = dec2bin('Text message looooooooool!!!!!1111\n',8);
bin = reshape(a.',[],1)' - '0';

initial = [1 bin]
N = length(initial);

%N = 1000000;
%initial = randi([0 1], N, 1)'; % генерация бит

bits = zeros(1,N);
bits(1) = initial(1);
for i=2:N
    bits(i)=mod(bits(i-1)+initial(i),2);
end

modData = mapping(bits, 1); % модуляция BPSK
nsignal = AWGnoise(modData, 20, 1);

offset = 1;
df = pi/5;
delta = complex(1, tan(df));
delta = delta / abs(delta);
for i=1:length(nsignal)
    nsignal(i) = nsignal(i) * offset;
    offset = offset * delta;
end

scatterplot(nsignal)

samples = nsignal;
demod = zeros(1,N);
demod(1)=1;
for i=2:N
    ang = mod(abs(angle(samples(i)) - angle(samples(i-1))), 2*pi);
    if (ang > (pi/2)) & (ang < (3*pi/2))
        demod(i) = 1;
    end
    
    if demod(i) ~= initial(i)
        i
        samples(i)
        samples(i-1)
        mod(ang, 2*pi)
        initial(i)
        mod(ang, 2*pi) > (pi/2)
    end
end

final = demod(2:end);

b = reshape(final, 8,[]).';
msg = num2str(b);
fprintf(char(bin2dec(msg)))

err = sum(xor(bin, final))
ber = err / N;