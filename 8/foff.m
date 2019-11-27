% Востановление несущей
%> @file foff.m
% =========================================================================
%> @brief Востановление несущей
%> @param fsymbols символы сигнала с частотным и фазовым сдвигом сдвигов с
%> с пилотами
%> @param step растояние между пилотными символами
%> @param pilots пилотные символы
%> @return symbols символы сигнала без частотного, фазового сдвига, пилотов
% =========================================================================
function symbols = foff (fsymbols, step, pilots)
% =========================================================================
    sk = [];
    leng = length(pilots);
    for i=1:leng+step:length(fsymbols)
        sk = [sk; fsymbols(i:i+leng-1)];
    end
            
    ck = sk.*conj(pilots);
    
    Rn(1:18) = 0;
    for i=1:18
        s = 0;
        
        for j=1:36-i
            for k = 1:300
                s = s + ck(k,j+i)*conj(ck(k,j));
            end
        end
        
        Rn(i)=s;
    end
    
    df = 2/(j*(j+1))*sum(angle(Rn))/(2*pi);
    
    fsymbols = fsymbols.*exp(sqrt(-1)*(-2*pi*[1:length(fsymbols)]*(df)));
    
    leng = length(pilots);
    symbols = [];
    for i=1:leng+step:length(fsymbols)
        symbols = [symbols fsymbols(i+leng:i+leng+step-1)];
    end

    leng = length(pilots);
    symbols = [];
    for i=1:leng+step:length(fsymbols)
        sk = [fsymbols(i:i+leng-1)];
        ck = sk.*conj(pilots);
        fi = angle(sum(ck));
        ex = exp(sqrt(-1)*(-fi));
        
        fsymbols = fsymbols*ex;
        symbols = [symbols fsymbols(i+leng:i+leng+step-1)];
    end
end
