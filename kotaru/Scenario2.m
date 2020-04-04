function Scenario2(runParam)
% The second (successful) attempt to create proper scenario for our purpose

%============== Создаем пакет ===========%
if runParam.subcarrierMode <= 3
    chanBW      = 'CBW40';            % Channel bandwidth MHz
else
    chanBW      = 'CBW20';            % Channel bandwidth MHz
end

mcs         = 0;                  % MCS: 0 - BPSK 1/2
psduLength  = 2048;
numTransmitAntennas = 1;

%Create a wlanHTConfig object for OFDM operation
cfgHT = wlanHTConfig('ChannelBandwidth', chanBW, ...
    'MCS', mcs, 'PSDULength', psduLength, 'NumTransmitAntennas', numTransmitAntennas);

fs = wlanSampleRate(cfgHT);  % Частота дискретизации

% Сгенерируем данные для передачи
psdu = randi([0 1], cfgHT.PSDULength*8, 1);

% Поле данных
HTData = wlanHTData(psdu,cfgHT);
size(HTData);

% Поля преамбулы
lstf = wlanLSTF(cfgHT);
lltf = wlanLLTF(cfgHT);
lsig = wlanLSIG(cfgHT);

% Объединяем поля в сообщение
txPPDU = [lstf; lltf; lsig; HTData];
sigLen = size(txPPDU, 1);
preamLen = size([lstf; lltf; lsig],1);

%============== Настраиваем канал ===========%

freq = runParam.fc;    %signal frequency
separation = 299792458/freq/2;  % lambda/2 - distance between antennas

% Настройка антенн
AA(1) = winner2.AntennaArray();  % BS, передатчик, 1 дефолтная антенна
AA(2) = winner2.AntennaArray('ULA', runParam.antennaNum, separation); % MS, приемник, 4 антенны

% theta = 0.5*pi();
% placement = [-1.5*separation*cos(theta) -1.5*separation*sin(theta) 0;
%     -0.5*separation*cos(theta) -0.5*separation*sin(theta) 0;
%     0.5*separation*cos(theta) 0.5*separation*sin(theta) 0;
%     1.5*separation*cos(theta) 1.5*separation*sin(theta) 0];
% 
% AA(2) = winner2.AntennaArray('Pos', placement);

arrayAngle = runParam.arrayAngle;
save('channelData.mat', 'arrayAngle');
AA(2).Rot = pi*[0 0 -arrayAngle/180].';

% Отрисовка массива не учитывает задаваемый поворот, т.к. это поворот всего
% массива как единого целого, а отрисовка это графический вывод координат в
% системе отсчета центра массива

% figure
% hold on
% pos = {AA(2).Element(:).Pos};
% plot(cellfun(@(x) x(1),pos),cellfun(@(x) x(2),pos),'+');
% xlim([-1 1]);
% ylim([-1 1]);
% title('AA Element Positions');
        
rndSeed = runParam.rndSeed;
rmax=50;    %Maximum layout range, specified as a scalar representing
            % the maximum layout range in meters used to randomly generate the MS and BS positions.


% Настройка сценария
% Информация: winner2.layoutparset(msIdx,bsIdx,K,arrays)
msIdx = 2;
bsIdx = {1};
numLinks = 1;
cfgLayout = winner2.layoutparset(msIdx, bsIdx, numLinks, AA, rmax, rndSeed);
cfgLayout.Pairing = [1; 2];
cfgLayout.ScenarioVector = 1;           % A1 scenario
cfgLayout.PropagConditionVector = runParam.isLos;    % 0 - NLOS, 1 - LOS

% Памятка:
% cfg.Layout.Stations(1) - это BS, передатчик, 1 антенна
% cfg.Layout.Stations(2) - это MS, приемник, 4 антенны

velocity = 1e-6;
cfgLayout.Stations(2).Velocity = [-velocity 0 0].';

% Можно самому выставить позиции, это будет влиять лишь на LOS angle.
cfgLayout.Stations(1).Pos= runParam.BSposition;   % BS position
cfgLayout.Stations(2).Pos = runParam.MSposition;   % MS position

% Set up model parameters for WINNER II channel
cfgModel = winner2.wimparset;
cfgModel.UseManualPropCondition = 'yes'; % If ‘yes’ the propagation condition (los/nlos) setting is defined manually in LAYOUTPARSET in PropagConditionVector.
cfgModel.CenterFrequency      = freq;   % Carrier frequency in Hz
cfgModel.RandomSeed         = rndSeed;    % Repeatability
cfgModel.UniformTimeSampling = 'yes';   % all links to be sampled at the same time instants

% Oversampling factor, number of time samples per half wavelength. For successful Doppler analysis, one should select SampleDensity > 1.
 cfgModel.SampleDensity = round(physconst('LightSpeed')/ ...
    freq/2/(velocity/fs));

if runParam.FixedPdpUsed
    cfgModel.FixedPdpUsed = 'yes'; % If ‘yes’ the power and delay parameters are not drawn randomly, but taken from the CDL parameter tables
else
    cfgModel.FixedPdpUsed = 'no';
end

if runParam.FixedAnglesUsed
    cfgModel.FixedAnglesUsed    = 'yes'; % If ‘yes’ the angle parameters are not drawn randomly, but taken from the CDL parameter tables
else
    cfgModel.FixedAnglesUsed    = 'no';
end

if runParam.IntraClusterDsUsed
    cfgModel.IntraClusterDsUsed = 'yes';    % If ‘yes’ the two strongest clusters in power are divided in delay into three subclusters.
else
    cfgModel.IntraClusterDsUsed = 'no';
end

if runParam.PolarisedArrays
    cfgModel.PolarisedArrays    = 'yes';    % Use dual-polarized arrays
else
    cfgModel.PolarisedArrays    = 'no';
end



% % По дефолту выключены, в примере для Wi-Fi не включались

if runParam.PathLossModelUsed
    cfgModel.PathLossModelUsed = 'yes';
else
    cfgModel.PathLossModelUsed = 'no';
end
    
if runParam.ShadowingModelUsed
    cfgModel.ShadowingModelUsed = 'yes'; 
else
    cfgModel.ShadowingModelUsed = 'no';
end


% Transmit through the WINNER II channel with 10 all-zero
% samples appended to account for channel filter delay
numPadZeros = 10;
cfgModel.NumTimeSamples = sigLen + numPadZeros; 

%============== Прогоняем через канал ===========%

save('runParam.mat', 'runParam');

% Create the WINNER II channel System object
WINNERChan = comm.WINNER2Channel(cfgModel, cfgLayout);

% Call the info method to check some derived channel parameters
%chanInfo = info(WINNERChan) %#ok<NOPTS>

if runParam.printChanInfo
    chanInfo = info(WINNERChan)
end

% Sound the WINNER II channel for all users
chanOut = WINNERChan([txPPDU;zeros(numPadZeros,1)]);

% % Определяем длину преамбулы в семплах
% fieldInd = wlanFieldIndices(cfgNHT);
% numPream = fieldInd.LSIG(2);
% 
% plotNum = 4;
% 
% % Выводим на график преамбулу, разделенную на L-STF, L-LTF, L-SIG
% figure
% hold on;
% time = ((0:double(numPream)-1)/fs)*1e6;   % Массив временных отметок в мкс
% peak = 1.2*max(abs(chanOut{1}(1:numPream, 1)));
% fieldMarkers = zeros(numPream,1);
% fieldMarkers(fieldInd.LSTF(2)-1,1)  = peak;
% fieldMarkers(fieldInd.LLTF(2)-1,1) = peak;
% fieldMarkers(fieldInd.LSIG(2)-1,1) = peak;
% plot(time,abs(chanOut{1}(1:numPream,1:plotNum)),time,fieldMarkers)
% xlabel ('Time (microseconds)')
% ylabel('Magnitude')
% title('Non-HT Format Preamble')

%============== Импульсная характеристика ===========%

% release(WINNERChan);
% 
% frameLen = cfgModel.NumTimeSamples;
% impulseSigTx = [ones(1,1);zeros(frameLen-1,1)];
% impulseSigRx = WINNERChan(impulseSigTx);
%
% % График импульсной характеристики
% figure
% hold on;
% stem((0:frameLen-1)/chanInfo.SampleRate, abs(impulseSigRx{1}(:,:)));
% minX = -0.1e-6;
% maxX = 2e-6;
% xlim([minX, maxX]);
% xlabel('Time (s)'); ylabel('Magnitude');
% title('Impulse Response for 4 antennas');
%     
%============== Обработка данных ===========%

% Наложение шума, в примере с мобильными станциями это не делается
% % Add AWGN
% rxSig = cellfun(@awgn, chanOut, ...
%     num2cell(snr*ones(numUsers, 1)), 'UniformOutput', false);

% mat = cell(numUsers,1);
% for uIdx = 1:numUsers
%     % Compute the feedback matrix based on received signal per user
%     mat{uIdx} = vhtCSIFeedback(rxSig{uIdx}(chanDelay(uIdx)+1:end,:), ...
%         cfgVHTNDP, uIdx, numSTSVec);
% end

signal = chanOut{1};
fieldInd = wlanFieldIndices(cfgHT);

save('signal.mat', 'signal');
save('signal.mat', 'fieldInd', '-append');

% from SpotFi
% SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26, ...
%     -22,-18,-14,-10, -6,-2,2,6,10,14,18,22,26,30,34, ...
%     38,42,46,50,54, 58].'; % WiFi subcarrier indices at which CSI is available

if runParam.subcarrierMode <= 3
    csi = zeros(114,runParam.antennaNum);
    for i=1:runParam.antennaNum
        antennaSignal = signal(:,i);
        doubleCSI = wlanHTLTFDemodulate(antennaSignal(fieldInd.LLTF(1):fieldInd.LLTF(2)), cfgHT);
        csi(:,i) = doubleCSI(:,1);
    end
    
    subCarr = [csi(1:57,:); zeros(3,runParam.antennaNum); csi(58:114,:)]; % zeros(3,4): -1, 0, 1, которых нет
else
    if runParam.doubleCSI
        csi = zeros(52,runParam.antennaNum,2);
        for i=1:runParam.antennaNum
            antennaSignal = signal(:,i);
            doubleCSI = wlanLLTFDemodulate(antennaSignal(fieldInd.LLTF(1):fieldInd.LLTF(2)), cfgHT);
            csi(:,i,1) = doubleCSI(:,1);
            csi(:,i,2) = doubleCSI(:,2);
        end
    else
        csi = zeros(52,runParam.antennaNum);
        for i=1:runParam.antennaNum
            antennaSignal = signal(:,i);
            doubleCSI = wlanLLTFDemodulate(antennaSignal(fieldInd.LLTF(1):fieldInd.LLTF(2)), cfgHT);
            csi(:,i) = doubleCSI(:,1);
        end
    end
    
    subCarr = [csi(1:26,:,:); zeros(1,runParam.antennaNum, size(csi,3)); csi(27:52,:,:)];
end

switch runParam.subcarrierMode
    case 0
        subCarrInd = (-58:4:58)+59;
    case 1
        subCarrInd = (-58:4:58)+59;
        subCarrInd = [subCarrInd(1:14) subCarrInd(17:end)]; % to remove -2, 2
    case 2
        subCarrInd = [-58:-2 2:58]+59;
    case 3
        subCarrInd = (-58:58)+59;
        subCarrInd = [subCarrInd(1:26) subCarrInd(28:53) subCarrInd(65:90) subCarrInd(92:end)];
    case 4
        subCarrInd = (-26:26)+27;
        subCarrInd = [subCarrInd(1:26) subCarrInd(28:end)];
end

rareCSI = subCarr(subCarrInd,:,:);

CSI = [];
for i=1:runParam.antennaNum
    CSI = [CSI; rareCSI(:,i,:)];
end

% old way just to check how it works
% vectorCSI = [];
% for i=1:runParam.packetNum
%     vectorCSI = [vectorCSI CSI];
% end

if runParam.doubleCSI
    vectorCSI = [CSI(:,:,1) CSI(:,:,2)];
else
    vectorCSI = CSI;
end

save('vectorCSI.mat', 'vectorCSI');

end