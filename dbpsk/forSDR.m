scatterplot(samples)

N = length(samples);

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

demod

final = demod(2:end);   % поменять пределы для вычленения сообщения из потока

b = reshape(final, 8,[]).';
msg = num2str(b);
fprintf(char(bin2dec(msg)))

err = sum(xor(bin, final))
ber = err / N;