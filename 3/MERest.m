% Вычисление MER 
%> @file MERest.m
% =========================================================================
%> @brief Вычисление MER
%> @param signal идеальный сигнал
%> @param Nsignal зашемленные сигнал
%> @return MER ошибка модуляции в дБ
% =========================================================================
function MER = MERest (signal, Nsignal)
    dif = Nsignal - signal;
    noise = sum(dif.*conj(dif))/length(dif);
    sig = sum(signal.*conj(signal))/length(signal);
    MER = 10 * log10(sig / noise);
end
