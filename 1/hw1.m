% Домашнее задание №1. Отображение бит на созвездие.
%> @file hw1.m
% Используемые файлы: mapping.m (рыба для функции отображения бит на созвездие).
% Задание: написать не используя стандартные функции и объекты, такие как 
% pskmod, apskmod, qammod,comm.[тип модуляции]Modulator и т.д.
% Используйте файл mapping.m, как основу для функции.
% Данный файл служит для тестирования вашей функции.
% Внимание: во избжание ошибок обратите внимание на разницу результатов
% функций dec2bin и de2bi. Считаем, что старший бит в передаваемой 
% посделовательности записан первым. 
% Для прохождения тестов нумерация созвездий должна совпадать с
% прописаннной в задании
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
%> Проверка созвездий
% =========================================================================
% Проверка осуществляется стандартными функциями демодуляции
for constellation = 1:5
    switch (constellation)
        case 1 % BPSK  
            BitInSym = 1;                 % колличество бит на точку
            ConstName = 'BPSK';
            demod = comm.BPSKDemodulator; % стандартный демодулятор
        case 2 % QPSK
            BitInSym = 2;                 % колличество бит на точку
            ConstName = 'QPSK';
            demod = comm.QPSKDemodulator; % стандартный демодулятор
        case 3 % 8PSK
            BitInSym = 3;                 % колличество бит на точку
            ConstName = '8PSK';
            demod = comm.PSKDemodulator('ModulationOrder', 8); % стандартный демодулятор
        case 4 % 16APSK
            BitInSym = 4;
            ConstName = '16APSK';
            % до версии 2018а нет стандартного демодулятора, поэтому тут
            % пока проверим визуально
        case 5 % 16QAM
            BitInSym = 4;
            ConstName = '16QAM';
    end
%> Проверяем расположение точек на созвездии
    data = (0:2^BitInSym-1);    % массив десятичных чисел
    bits = de2bi(data, BitInSym);   % матрица, строки которой являются
    % бинарным представлением вышеупомянутых чисел

    % .' - транспонирование
    % reshape - изменение размеров матрицы с сохранением элементов
    % reshape(..., 1, []) - хотим матрицу с 1 строкой и не паримся о
    % количестве столбцов
    
    %bits = reshape(bits(:,end:-1:1).', 1, []);
    bits = bits(:,end:-1:1);
    bits = bits.';
    bits = reshape(bits, 1, []);
    % в итоге получили поток битов, который состоит из наших бинарных чисел
    
    % кидаем этот поток в нашу функцию, которая, видимо, возвращает
    % массив комплексных чисел
    modData = mapping(bits, constellation);
    
    %> Визуализация
    scatterplot(modData)
    text(real(modData)+0.1, imag(modData), dec2bin(data))
    title(ConstName)
    axis([-2 2 -2 2])
    %> Считаем среднюю мощность (проверка нормировки)
    P = sum(modData.*conj(modData))/length(modData);
    if abs(1-P)>0.00001
        Error = 'Проверьте нормировку созвездия'
        ConstName = ConstName
    end
    if constellation < 4
        %> Демодулируем стандартным демодулятором
        bits = randi([0 1], 1, 120000); % генерация бит
        modData = mapping(bits, constellation);
        checkData = demod(modData.');
        checkBits = de2bi(checkData, BitInSym);
        checkBits = reshape(checkBits(:,end:-1:1).', 1, []);
        Nerr = sum(xor(checkBits,bits))
        if Nerr~= 0
            Error = 'Проверьте созвездие'
            ConstName = ConstName
        end
    elseif constellation == 5
        %> Демодулируем стандартным демодулятором
        bits = randi([0 1], 1, 120000); % генерация бит
        modData = mapping(bits, constellation);
        norm = sqrt(10);
        modData = norm*modData;
        %modData = qammod(bits, 16);
        checkData = qamdemod(modData.', 16);
        %checkData = qamdemod(modData, 16);%.', 16);
        checkBits = de2bi(checkData, BitInSym);
        checkBits = reshape(checkBits(:,end:-1:1).', 1, []);
        Nerr = sum(xor(checkBits,bits))
        if Nerr~= 0
            Error = 'Проверьте созвездие'
            ConstName = ConstName
        end
    end
end

    
