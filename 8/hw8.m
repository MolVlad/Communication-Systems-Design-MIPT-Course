% Домашнее задание №8. Частотная синхргонизация.
%> @file hw8.m
%> Используйте функции: mapping.m, AWGnoise.m
% =========================================================================
%> Подготовка рабочего места
% =========================================================================
    %> Отчистка workspace
    clear all;
    %> Закрытие рисунков
    close all;
    %> Отчистка Command Window
    clc;
% =========================================================================
% Задание: Напишите функцию частотной синхронизации
% =========================================================================
N = 2 * 90 * 16 * 300;
bits = randi([0 1], N, 1); % генерация бит
modData = mapping(bits, 2); % модуляция
% вставляем пилоты
pilots(1:36) = exp(sqrt(-1)*pi/4);
step = 90*16;
leng = length(modData);

pmodData = [];
for i = 1:step:leng
    pmodData = [pmodData, pilots, modData(i:i+step-1)];
end
SNR = 20;
nsignal = AWGnoise (pmodData, SNR, 1); % функция, написанная вами в hw2
scatterplot(nsignal)
freqoff = 0.0005%*rand(1) % нормализованный частотный сдвиг oт 0 до 1
phaseoff = 2*pi*rand(1);
%> Частотный сдвиг <=> поворот
leng = length(nsignal);
nsignal = nsignal.*exp(sqrt(-1)*(2*pi*[1:leng]*freqoff+phaseoff));
symbols = foff (nsignal, step, pilots);
scatterplot(symbols)
mer = comm.MER;
MER = mer(modData.', symbols.')
if abs (MER-SNR)>3
    'Проверьте функцию'
end
