% Генерация коэффицентов (импульсной характеристики)для фильтра 
% приподнятого косинуса
%> @file RCcoeff.m
% =========================================================================
%> @brief Генерация коэффицентов (импульсной характеристики)для фильтра 
%> приподнятого косинуса
%> @param span Длина фильтра в символах (число боковых лепестков sinc, сумма
%> с двух сторон)
%> @param nsamp Число выборок на символ
%> @param rolloff Коэффицент сглаживания (alfa)
%> @return coeff коэффиценты для фильтра корень из приподнятого косинуса
% =========================================================================
function coeff = RCcoeff (span, nsamp, rolloff)
    coeff = zeros(1, nsamp * span + 1);
    
    T = nsamp;
    t0 = span /2 * nsamp + 1;  % момент времени t=0
    t1 = T/(2 * rolloff);
    
    coeff(t0) = 1/sqrt(T);
    
    for i=1:(span * nsamp / 2)
        if i ~= t1
            coeff(t0 + i) = coeff(t0) * sin(pi * i / T) * cos(pi * rolloff * i / T);
            coeff(t0 + i) = coeff(t0 + i) / (pi * i / T) / (1 - (2 * rolloff * i / T) ^ 2);
            coeff(t0 - i) = coeff(t0 + i);
        else
            coeff(t0 + t1) = coeff(t0) * rolloff / 2 * sin(pi/(2 * rolloff));
            coeff(t0 - t1) = coeff(t0 + t1);            
        end
    end  
end
