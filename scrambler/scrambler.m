%> Отчистка workspace
clear all;
%> Закрытие рисунков
close all;
%> Отчистка Command Window
clc;

initSeq = [1 0 0 1 0 1 0 1 0 0 0 0 0 0 0];

%==========================================
% Transmitter
%==========================================

%txBits = randi([0 1], 1, 1200); % генерация бит
txBits = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

txSeq = initSeq;

txData(1:length(txBits)) = 0;
for i=1:length(txBits)
    res = xor(txSeq(end),txSeq(end-1));
    txData(i) = xor(txBits(i),res);
    txSeq = [res txSeq(1:end-1)];
end

txBits
txData

modData = mapping(txData, 1); % модуляция

%==========================================
% Receiver
%==========================================

rxData = demapping(modData, 1, 0)';

rxSeq = initSeq;

rxBits(1:length(rxData)) = 0;
for i=1:length(rxData)
    res = xor(rxSeq(end),rxSeq(end-1));
    rxBits(i) = xor(rxData(i),res);
    rxSeq = [res rxSeq(1:end-1)];
end

ERR = Nerr(txBits, rxBits)
