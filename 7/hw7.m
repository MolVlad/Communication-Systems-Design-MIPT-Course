% Домашнее задание №7. Символьная синхронизация
%> @file hw7.m
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
% Задание 1. Написать генератор m-последовательности 63 (6, 5). 
% Начальное состояние буфера {1, 0, 0, 0, 0, 0}. Изучить ее
% авторорреляционные своства (модулировать BPSK, собрать кадр из 5 BPSK 
% последовательностей, посчитать корреляцию с исходной посдледовательностью
% по всей длине кадра, изобразить результат)
% =========================================================================
    seq = Mseq();
    
    n = length(seq);
    data = [seq seq seq seq seq];
    
    modSeq = mapping(seq,1);
    modData = mapping(data,1);
    
    r = [];
    for j=1:4*n+1
        sum=0;
        
        for i=1:n
            sum = sum + modData(i+j-1) * modSeq(i);
        end
        
        r = [r sum/63];
    end
    
    x = zeros(1, length(r));
    for i=1:length(r)
        x(i)=i;
    end

    f = figure;
    stem(x,abs(r));
    grid on;
    title('M-seq autocorrelation');
    name = 'autocorrelation.png';
    saveas(f, name);
    
% =========================================================================
% Задание 2. Собрать кадр seq + 900 случайных бит. 
% Собрать последовательность из 5 кадров. Наложить щум.
% Исследовать корреляционные свойства с исходной последовательностью
% от -5 до 20 дБ. Нарисовать графики. Сделать вывод о пороговом значении
% отношения сигнал-шум, при котором можно обнаружить пик.
% =========================================================================

frame = [randi([0 1], 1, 400) seq];
data = [frame frame frame frame frame];
modData = mapping(data, 1);

%SNR = -5:1:20;
SNR = [-2 5];

for k = 1:length(SNR)
    noiseData = AWGnoise(modData, SNR(k),1);
    
    r = [];
    for i=1:5*length(frame)-n+1
        sum=0;
        
        for j=1:n
            sum = sum + noiseData(i+j-1) * modSeq(j);
        end
        
        r = [r sum/63];
    end
    
    x = zeros(1, length(r));
    for i=1:length(r)
        x(i)=i;
    end

    f = figure;
    stem(x,abs(r));
    grid on;
    
    if SNR(k) == -2
        title('Signal correlation SNR = -2');
        name = 'signal_correlation_snr(-2).png';
    elseif SNR(k) == 5
    	title('Signal correlation SNR = 5');
        name = 'signal_correlation_snr(5).png';
    end
    
	saveas(f, name);

    % Пороговое значение порядка -2dB
end


% =========================================================================
% Задание 3. Собрать кадр seq + 900 случайных бит. 
% Собрать последовательность из 5 кадров. Наложить частотный сдвиг.
% Исследовать корреляционные свойства при частотном сдвиге от -0,2 до 0,2 pi.
% Частотный сдвиг реализовать по формуле из лекции. 
% Нарисовать графики. Сделать вывод о пороговом значении
% частотного сдвига, при котором можно обнаружить пик.
% =========================================================================

frame = [randi([0 1], 1, 900) seq];
data = [frame frame frame frame frame];
modData = mapping(data, 1);

rotateData(1:length(modData)) = 0;

freqOffset = [0.2 0.02];
freqOffset = freqOffset * pi;

for k = 1:length(freqOffset)
    
    offset = 1;
    delta = complex(1, tan(freqOffset(k)));
    delta = delta / abs(delta);
    
    for i=1:length(modData)
    	rotateData(i) = modData(i) * offset;
        offset = offset * delta;
    end
    
    r = [];
    for i=1:5*length(frame)-n+1
        sum=0;
        
        for j=1:n
            sum = sum + rotateData(i+j-1) * modSeq(j);
        end
        
        r = [r sum/63];
    end
    
    x = zeros(1, length(r));
    for i=1:length(r)
        x(i)=i;
    end

    f = figure;
    stem(x,abs(r));
	grid on;
     
    if freqOffset(k) == 0.2*pi
        title('Signal correlation freq offset = 0.2pi');
        name = 'signal_correlation_freqOffset(0.2pi).png';
    elseif freqOffset(k) == 0.02*pi
        title('Signal correlation freq offset = 0.02pi');
        name = 'signal_correlation_freqOffset(0.02pi).png';
    end
    
	saveas(f, name);

    %Пороговое значение порядка 0.02pi
end

% =========================================================================
% Задание 4. Проделать 2 и 3 для дифференциальных коэффициентов.
% Коррелировать с дифф. коэф. от исходной последовательности. 
% =========================================================================

difModSeq(1, n) = 0;

for i=1:n-1
    difModSeq(i) = modSeq(i)*conj(modSeq(i+1));
    difModSeq(i) = difModSeq(i)/sqrt(modSeq(i)*conj(modSeq(i)*modSeq(i+1)*conj(modSeq(i+1))));
end

difModSeq(n) = modSeq(n)*conj(modSeq(1));
difModSeq(n) = difModSeq(n)/sqrt(modSeq(n)*conj(modSeq(n)*modSeq(1)*conj(modSeq(1))));

frame = [randi([0 1], 1, 400) seq];
data = [frame frame frame frame frame];
modData = mapping(data, 1);

%SNR = -5:1:20;
SNR = [0 5];

for k = 1:length(SNR)
    noiseData = AWGnoise(modData, SNR(k),1);
    
    l = length(noiseData);
    difNoiseData(1, l) = 0;

    for i=2:l-1
        difNoiseData(i) = noiseData(i)*conj(noiseData(i+1));
        difNoiseData(i) = difNoiseData(i)/sqrt(noiseData(i)*conj(noiseData(i)));
        difNoiseData(i) = difNoiseData(i)/sqrt(noiseData(i+1)*conj(noiseData(i+1)));
    end
    
    r = [];
    for i=1:5*length(frame)-n+1
        sum=0;
        
        for j=1:n
            sum = sum + difNoiseData(i+j-1) * difModSeq(j);
        end
        
        r = [r sum/63];
    end
    
    x = zeros(1, length(r));
    for i=1:length(r)
        x(i)=i;
    end

    f = figure;
    stem(x,abs(r));
    grid on;
    
    if SNR(k) == 0
        title('Dif signal correlation SNR = 0');
        name = 'dif_signal_correlation_snr(0).png';
    elseif SNR(k) == 5
    	title('Dif signal correlation SNR = 5');
        name = 'dif_signal_correlation_snr(5).png';
    end
    
	saveas(f, name);

    % Пороговое значение порядка 0dB
end

frame = [randi([0 1], 1, 900) seq];
data = [frame frame frame frame frame];
modData = mapping(data, 1);

rotateData(1:length(modData)) = 0;

freqOffset = [5 2];
freqOffset = freqOffset * pi;

for k = 1:length(freqOffset)
    
    offset = 1;
    delta = complex(1, tan(freqOffset(k)));
    delta = delta / abs(delta);
    
    for i=1:length(modData)
    	rotateData(i) = modData(i) * offset;
        offset = offset * delta;
    end
    
    l = length(rotateData);
    difRotateData(1, l) = 0;

    for i=2:l-1
        difRotateData(i) = rotateData(i)*conj(rotateData(i+1));
        difRotateData(i) = difRotateData(i)/sqrt(rotateData(i)*conj(rotateData(i)));
        difRotateData(i) = difRotateData(i)/sqrt(rotateData(i+1)*conj(rotateData(i+1)));
    end
    
    r = [];
    for i=1:5*length(frame)-n+1
        sum=0;
        
        for j=1:n
            sum = sum + difRotateData(i+j-1) * difModSeq(j);
        end
        
        r = [r sum/63];
    end
    
    x = zeros(1, length(r));
    for i=1:length(r)
        x(i)=i;
    end

    f = figure;
    stem(x,abs(r));
	grid on;
     
    if freqOffset(k) == 5*pi
        title('Dif signal correlation freq offset = 5pi');
        name = 'dif_signal_correlation_freqOffset(5pi).png';
    elseif freqOffset(k) == 2*pi
        title('Dif signal correlation freq offset = 2pi');
        name = 'dif_signal_correlation_freqOffset(2pi).png';
    end
    
	saveas(f, name);

    %Пороговое значение без шума может быть сколь угодно большим
end

% =========================================================================
% Для проверки: прислать
% 1. Эти файлы
% 2. Картинку с автокорреляцией
% 3. По 2 картинки с шумами (указать значение) для обычной корреляции и для дифф.
% 4. По 2 картинки с частотными сдвигами (указать значение) для обычной корреляции и для дифф.
% 5. Выводы о пороговых значениях.
% =========================================================================
