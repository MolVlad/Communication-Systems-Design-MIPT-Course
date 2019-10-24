% Домашнее задание №4. Канальное кодирование. Кодирование и декодирование
%> @file hw4.m
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
%> Задание 1. Реализуйте алгоритм итеративного декодирование для примера из лекции.
%> Смысл - понять что такое пересылка сообщений.
% =========================================================================
%> Принятое слово
% C = [1 1 0 1 0 1 0 1];
% %> Проверочная матрица
% H = [0 1 0 1 1 0 0 1;...
%      1 1 1 0 0 1 0 0;...
%      0 0 1 0 0 1 1 1;...
%      1 0 0 1 1 0 1 0];
% %> bits - переданное сообщение. Нужно вычислить
% %> @todo место для вашго кода
% 
f = zeros(1,4); % check nodes
f(1) = 1;   % for start while loop
bits1 = C;     %variable nodes
attempts = 0;
threshold = 5;

while (sum(f) > 0 & attempts < threshold)
    attempts = attempts + 1;
    f = zeros(1,4);
    for i=1:4
        for j=1:8
            if H(i,j) == 1
                f(i) = f(i) + bits1(j);
            end
        end
    end
    
    f = mod(f,2)

    buf(1:8) = C;
    
    for j=1:8
        for i=1:4
            if H(i,j) == 1
                if f(i) == 0
                    buf(j) = buf(j) + bits1(j);
                else
                    buf(j) = buf(j) + mod(bits1(j) + 1, 2);  % invert
                end
            end
        end
        if buf(j) > 1
            bits1(j) = 1;
        else
            bits1(j) = 0;
        end
    end    
end

if sum(f) > 0
    'unsuccessful'
end

% =========================================================================
chack = [1 0 0 1 0 1 0 1];
err = sum(xor(chack, bits1))

% =========================================================================
%> Задание 3. 
%> Использюя стандартный кодер и декодер LDPC кода comm.LDPCEncoder,
%> comm.LDPCDecoder и ваши функции mapping demapping для QPSK постройте BER 
%> характеристку до 10^-7 для кода dvbs2ldpc(1/2) (default) для 10 и 50 итераций
%> декодирования. Напищите в комментариях какое решение дает лучшй
%> результат, почему? 
% =========================================================================
%> @todo место для вашго кода

LDPCenc = comm.LDPCEncoder; %> инициализация кодера
LDPCdec50 = comm.LDPCDecoder; %> инициализация декодера 50 итераций (дефолт)
LDPCdec10 = comm.LDPCDecoder; %> инициализация декодера для 10 итераций
LDPCdec10.MaximumIterationCount = 10;   % 10 итераций

constellation = 2;

N = 32400;
BER50 = [];
SNR50 = [0:0.05:1];

for i=1:length(SNR50)
    Ne50 = 0;
    for frame = 1:50 % номер блока кода 
            i
            frame
            
            bits50 = randi([0 1], N, 1); % генерация бит
            LDPCcod50 = LDPCenc(bits50); % кодирование
            modData50 = mapping(LDPCcod50.', constellation); % модуляция
            nsignal50 = AWGnoise (modData50, SNR50(i), 1); % функция, написанная вами в hw2
            LLR50 = demapping (nsignal50, constellation, 1, SNR50(i)); % мягкое решение\
            checkBits50 =  LDPCdec50(LLR50);
            ERR50 = Nerr (bits50, checkBits50);
            Ne50 = Ne50+ERR50;
    end
    
    BER50 = [BER50, Ne50/(N*frame)];
end

BER10 = [];
SNR10 = [0:0.05:3];

for i=1:length(SNR10)
    Ne10 = 0;
    for frame = 1:50 % номер блока кода 
            i
            frame
            
            bits10 = randi([0 1], N, 1); % генерация бит            
            LDPCcod10 = LDPCenc(bits10); % кодирование
            modData10 = mapping(LDPCcod10.', constellation); % модуляция
            nsignal10 = AWGnoise (modData10, SNR10(i), 1); % функция, написанная вами в hw2
            LLR10 = demapping (nsignal10, constellation, 1, SNR10(i)); % мягкое решение\
            checkBits10 =  LDPCdec10(LLR10);
            ERR10 = Nerr (bits10, checkBits10);
            Ne10 = Ne10+ERR10;
    end
    
    BER10 = [BER10, Ne10/(N*frame)];
end

f = figure;
semilogy(SNR50, BER50, 'r-')
hold on
semilogy(SNR10, BER10, 'b-')
grid on
title('BPSK LDPC two number of iterations')
xlabel('SNR(dB)')
ylabel('BER')
legend('50 iteration','10 iterations')
name = 'twoIterationNumber.png';
saveas(f, name);
name2 = 'twoIterationNumber.fig';
savefig(f,name2);

% Вывод:
% Декодер с максимальным числом итераций 50 дал более хороший результат,
% чем декодер с максимальным числом итераций 10, потому что за бОльшее
% число итераций он мог исправлять бОльшее число ошибок и сохранять низкий
% BER даже при низких значениях SNR.
% Однако он намного более затратный в плане вычислительных мощностей

% =========================================================================
%> Задание 4. 
%> Использюя стандартный кодер и декодер LDPC кода comm.LDPCEncoder,
%> comm.LDPCDecoder и ваши функции mapping demapping постройте BER 
%> характеристку до 10^-7 для двух кодов из dvbs2ldpc().
%> В комментариях расчитайте предел Шеннона для данной кодовой скорости и 
%> типа модуляции,сравните со значением по уровню 10^-7 (оцененым из графика). 
%> Оцените качество кода.
% =========================================================================
%> @todo место для вашго кода

constellation = 2;

p1 = dvbs2ldpc(9/10);
LDPCenc1 = comm.LDPCEncoder(p1); %> инициализация кодера 1
LDPCdec1 = comm.LDPCDecoder(p1); %> инициализация декодера 1

N1 = size(p1,2) - size(p1,1)
BER1 = [];
SNR1 = [-5:0.05:7];

for i=1:length(SNR1)
    Ne1 = 0;
    for frame = 1:50 % номер блока кода 
            i
            frame
            
            bits1 = randi([0 1], N1, 1); % генерация бит
            LDPCcod1 = LDPCenc1(bits1); % кодирование
            modData1 = mapping(LDPCcod1.', constellation); % модуляция
            nsignal1 = AWGnoise (modData1, SNR1(i), 1); % функция, написанная вами в hw2
            LLR1 = demapping (nsignal1, constellation, 1, SNR1(i)); % мягкое решение\
            checkBits1 =  LDPCdec1(LLR1);
            ERR1 = Nerr (bits1, checkBits1);
            Ne1 = Ne1+ERR1;
    end
    
    BER1 = [BER1, Ne1/(N1*frame)];
end

p2 = dvbs2ldpc(1/4);
LDPCenc2 = comm.LDPCEncoder(p2); %> инициализация кодера 2
LDPCdec2 = comm.LDPCDecoder(p2); %> инициализация декодера 2

N2 = size(p2,2) - size(p2,1)
BER2 = [];
SNR2 = [-5:0.05:-2.4];

for i=1:length(SNR2)
    Ne2 = 0;
    for frame = 1:50 % номер блока кода 
            i
            frame
            
            bits2 = randi([0 1], N2, 1); % генерация бит            
            LDPCcod2 = LDPCenc2(bits2); % кодирование
            modData2 = mapping(LDPCcod2.', constellation); % модуляция
            nsignal2 = AWGnoise (modData2, SNR2(i), 1); % функция, написанная вами в hw2
            LLR2 = demapping (nsignal2, constellation, 1, SNR2(i)); % мягкое решение\
            checkBits2 =  LDPCdec2(LLR2);
            ERR2 = Nerr (bits2, checkBits2);
            Ne2 = Ne2+ERR2;
    end
    
    BER2 = [BER2, Ne2/(N2*frame)];
end

f = figure;
semilogy(SNR1, BER1, 'r-')
hold on
semilogy(SNR2, BER2, 'b-')
grid on
title('BPSK LDPC two code speeds')
xlabel('SNR(dB)')
ylabel('BER')
legend('dvbs2ldpc(9/10)','dvbs2ldpc(1/4)')
name = 'twoCodesCompare.png';
saveas(f, name);
name2 = 'twoCodesCompare.fig';
savefig(f,name2);

% Вывод:

% Предел Шеннона: C/B=log2(1+SNR)
% где C - предельная пропускная способность канала (бит/c),
% B - полоса сигнала (Гц), SNR - отношение сигнал/шум.
% SNR = Psignal/Pnoise = (Eb*C)(No*B), (*) где
% Eb - энергия на 1 бит, No - спектральная плотность шума.
% Тогда предел Шеннона можно записать как:
% C/B=log2(1+Eb/No*(C/B)), откуда Eb/No=B/C*(2^(C/B)-1).

% Заметим, что C выражается как произведение символьной скорости,
% скорости кода и количества бит в символе.
% Символьная скорость выражается как 1/T, где T - период следования
% символов, причем для идеального формирующего фильтра 1/T=B, значит
% соотношение B/C просто равняется произведению скорости кода на количество
% бит на символ.

% В нашей лабораторной работе мы использовали BPSK сигналы, для которых
% количество бит на символ равняется 1. Таким образом, используя значение
% скорости кода мы можем рассчитать предел Шеннона для используемой
% сигнально-кодовой конструкции (код+модуляция).

% Eb/No - отношение энергии на один бит к спектральной мощности шума.
% Для кода dvbs2ldpc(1/4) предел Шеннона:
% Eb/No = 10*log10(4/1*(2^(1/4*1)-1))= -1.2100 dB.
% Для кода dvbs2ldpc(9/10):
% Eb/No = 10*log10(10/9*(2^(9/10*1)-1))=-0.1669 dB. 

% Теперь из полученных графиков зависимости BER от SNR мы можем определить,
% насколько наша сигнально-кодовая конструкция близка к теоретическому
% пределу. Для этого зададим приемлемый уровень ошибок - BERmin=10^(-5).
% Из графика определяем, какой минимальный уровень SNR нужно обеспечить,
% чтобы BER был меньше этого значения. Затем используя формулу (*)
% высчитываем Eb/No, соответствующий этому условию. 

% Для скорости кода 1/4 из графика SNR = -2.60419 dB,
% C/B=10*log10(1/4)=-6.0206 dB, Eb/No=-2.60419 - (-6.0206)=3.4164 dB.
% Для скорости кода 9/10 из графика SNR = 6.2839dB,
% C/B=10*log10(9/10)=-0.4576 dB, Eb/No=6.2839-(-0.4576)= 6.7415 dB.

% Итого:
% Для кода dvbs2ldpc(1/4) из полученных данных BER=10^(-5) при Eb/No= 3.4164 dB,
% а предел Шеннона для данной скорости кода (9/10): Eb/No= -1.2100 dB.
% Для кода dvbs2ldpc(9/10) из полученных данных BER=10^(-5) при Eb/No= 6.7415 dB,
% а предел Шеннона для данной скорости кода (9/10): Eb/No= -0.1669.

% Таким образом, расстояние до предела Шеннона составляет порядка 4-7 dB.