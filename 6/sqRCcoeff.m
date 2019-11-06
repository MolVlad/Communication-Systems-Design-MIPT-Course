% Генерация коэффицентов (импульсной характеристики)для фильтра корень 
% из приподнятого косинуса
%> @file sqRCcoeff.m
% =========================================================================
%> @brief Генерация коэффицентов (импульсной характеристики)для фильтра корень 
%> из приподнятого косинуса
%> @param span Длина фильтра в символах (число боковых лепестков sinc, сумма
%> с двух сторон)
%> @param nsamp Число выборок на символ
%> @param rolloff Коэффицент сглаживания (alfa)
%> @return coeff коэффиценты для фильтра корень из приподнятого косинуса
% =========================================================================
function coeff = sqRCcoeff (span, nsamp, rolloff)
    coeff = zeros(1, nsamp * span + 1);
    
    T = nsamp;
    t0 = span * nsamp / 2 + 1;  % момент времени t=0        
    t1 = T/(4 * rolloff);

    coeff(t0) = 1/sqrt(T) * ((1 - rolloff) + 4 * rolloff / pi);
    
    for i=1:(span * nsamp / 2)
        if i ~= t1
            coeff(t0 + i) = 1/sqrt(T)*(sin(pi*i*(1-rolloff)/T)+4*rolloff*i/T*cos(pi*i*(1+rolloff)/T));
            coeff(t0 + i) = coeff(t0 + i) / (pi * i / T) / (1 - (4 * rolloff * i / T) ^ 2);
            coeff(t0 - i) = 1/sqrt(T)*(sin(pi*(-1)*i*(1-rolloff)/T)+4*rolloff*(-1)*i/T*cos(pi*(-1)*i*(1+rolloff)/T));
            coeff(t0 - i) = coeff(t0 - i) / (pi * (-1) * i / T) / (1 - (4 * rolloff * (-1) * i / T) ^ 2);
        else
            coeff(t0 + t1) = rolloff/sqrt(2*T)*((1+2/pi)*sin(pi/(4*rolloff))+(1-2/pi)*cos(pi/(4*rolloff)));
            coeff(t0 - t1) = coeff(t0 + t1);
        end
    end
end
