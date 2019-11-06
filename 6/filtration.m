% Фильтрация
%> @file filtration.m
% =========================================================================
%> @brief Фильтрация
%> @param sign входной сигнал сигнал
%> @param coeff коэффиценты фильтра
%> @param nsamp число выборок на символ
%> @param UpSempFlag [1] -  фильтр с передескретезацией,[0] - фильтр без передескретизации 
%> @return filtsign отфильтрованный сигнал 
% =========================================================================
function filtsign = filtration(sign, coeff, nsamp, UpSempFlag)
    if UpSempFlag == 0
        nsamp = 1;
    end
    
    filtsign = zeros(1, nsamp * length(sign));
    
    span = (length(coeff) - 1) / 2;
    
    l = length(coeff);
    buf = zeros(1, l);
    
    for i = 1:length(sign)
         % умножаем, складываем в буфер
        buf = buf + coeff * sign(i);
        
        % записываем из буфера в filtsign
        filtsign(i*nsamp-nsamp+1:i*nsamp) = buf(1:nsamp);
        
        % сдвигаем буфер на nsamp
        buf = [buf(nsamp+1:end) zeros(1,nsamp)];
    end
end 
