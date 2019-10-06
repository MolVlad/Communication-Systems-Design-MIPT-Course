% Создаем конфигурационный объект для nonHT
nonHT = wlanNonHTConfig;

nonHT.NumTransmitAntennas = 4;
nonHT.MCS = 3;

% Конф объект для nonHT только с другой модуляцией
% DSSS - direct-sequence spread spectrum
nonHT2 = wlanNonHTConfig('Modulation','DSSS');

% Сгенерируем waveform для передачи
bits = [1;0;0;1];
txWaveform = wlanWaveformGenerator(bits,nonHT);