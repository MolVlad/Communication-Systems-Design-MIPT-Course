% Домашнее задание №3. Обратное отображение сигнала на созвездие, жесткое и мягкое решение. 
%> @file hw3.m
% Используемые файлы: mapping.m (написанная вами функция отображения бит на созвездие)
%                     AWGnoise.m(написанная вами функция добавления шума)
%                     Nerr.m    (написанная вами функция добавления шума)
% Задание:
% Написать функцию обратного отображения бит на созвездие demapping.m
% Расчитать водопадные характеристики для жестких решений
% до 10^(-5) для всех созвездий, сравнить с эталонными. 
% Расчитать водопадную характеристику до 10^(-5), для мягкого решения при
% LDPC кодировании. Сравнить с характеристикой без кодирования для QPSK
%> @warning тесты 1.2 2.2 - требуют большого колличества вычислительных
%> ресурсов (не обязательны к выполнению)
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
%> Задание 1. Жесткое решение
% =========================================================================
for constellation = 1:5
    % Выбор созвездия
    switch (constellation)
        case 1 % BPSK  
            BitInSym = 1;                 % колличество бит на точку
            ConstName = 'BPSK';
            SNR = (1:11);                 % Es/No дБ
        case 2 % QPSK
            BitInSym = 2;                 % колличество бит на точку
            ConstName = 'QPSK';
            SNR = (4:14); 
        case 3 % 8PSK
            BitInSym = 3;                 % колличество бит на точку
            ConstName = '8PSK';
            SNR = (8:18); 
        case 4 % 16APSK
            BitInSym = 4;
            ConstName = '16APSK';
            SNR = (11:21); 
        case 5 % 16QAM
            BitInSym = 4;
            ConstName = '16QAM';
            SNR = (11:21); 
    end
    % =====================================================================
    %> Проверка 1.1 Без шума
    % =====================================================================
    bits = randi([0 1], 1, 120000); % генерация бит
    modData = mapping(bits, constellation); % модуляция
    %> Функция которую вы должы написать
    checkBits = demapping (modData, constellation, 0);
    ERR = Nerr (bits, checkBits)
    if ERR == 0
        ans = 'Проверка 1.1 - успех'
        ConstName = ConstName
    else
        ans = 'Проверка 1.1 - провал'
        ConstName = ConstName
    end
%     % =====================================================================
%     %> Проверка 1.2 Водопадные характеристики
%     % =====================================================================
%     % Внимание: может считаться довольно долго!!!
%     BER = [];
%     for i = 1:length(SNR)
%         N = 3000000; % Число бит, если считет слишком долго - умешите.
%         bits = randi([0 1], 1, N); % генерация бит
%         modData = mapping(bits, constellation); % модуляция
%         nsignal = AWGnoise (modData, SNR(i), 1); % функция, написанная вами в hw2
%         chackBits = demapping (nsignal, constellation, 0);
%         ERR = Nerr (bits, chackBits)
%         BER = [BER, ERR/N];
%     end
%     f = figure;
%     semilogy(SNR, BER, 'r^-')
%     hold on
%     semilogy(SNR-10*log10(BitInSym), BER, 'bo-')
%     grid on
%     title(ConstName)
%     xlabel('EbNo(dB)/EsNo(dB)')
%     ylabel('BER')
%     %> Сохраните полученные изображения
%     ConstName = [ConstName,'.png'];
%     saveas(f, ConstName);
end
% =========================================================================
%> Проверка 2.1 Мягкое решение
% =========================================================================
    %> Проверка функции с мягким решением
    %> Стандартный LDPC код для длинного кадра (итоговая длина 64800 бит)
    %> стандарта DVB-S2 кодовая скорость 1/2, максимальное колличесво
    %> итераций 50
    LDPCenc = comm.LDPCEncoder; %> инициализация кодера
    LDPCdec = comm.LDPCDecoder; %> инициализация декодера
    N = 32400; % Число бит (по размеру матрицы LDPC кода)
    bits = randi([0 1], N, 1); % генерация бит
    LDPCcod = LDPCenc(bits); % кодирование
    modData = mapping(LDPCcod.', 5); % модуляция
    nsignal = AWGnoise (modData, 7, 1); % функция, написанная вами в hw2
    LLR = demapping (nsignal, 5, 1, 7); % мягкое решение
    checkBits =  LDPCdec(LLR);
    ERR = Nerr (bits, checkBits)
    if ERR == 0
        ans = 'Проверка 2.1 - успех'
    else
        ans = 'Проверка 2.1 - провал'
    end
    % =====================================================================
    %> Проверка 2.2 - Мягкое решение
    % =====================================================================
    BER = [];
    SNR = [5:0.05:7];
    for i = 1:length(SNR) % диапозон изменения SNR
        Ne = 0;
        for frame = 1:300 % номер блока кода 
            bits = randi([0 1], N, 1); % генерация бит
            LDPCcod = LDPCenc(bits); % кодирование
            modData = mapping(LDPCcod.', 5); % модуляция
            nsignal = AWGnoise (modData, SNR(i), 1); % функция, написанная вами в hw2
            LLR = demapping (nsignal, 5, 1, SNR(i)); % мягкое решение
            chackBits =  LDPCdec(LLR);
            ERR = Nerr (bits, chackBits);
            Ne = Ne+ERR
            frame
        end
        BER = [BER, Ne/(N*frame)];
    end
        f = figure;
    semilogy(SNR, BER, 'r^-')
    hold on
    semilogy(SNR-10*log10(2/2), BER, 'bo-')
    grid on
    title('QPSK LDPC 1/2')
    xlabel('EbNo(dB)/EsNo(dB)')
    ylabel('BER')
    %> Сохраните полученные изображения
    ConstName = ['QPSK_LDPC.fig'];
    savefig(f,ConstName);
    
    
    
