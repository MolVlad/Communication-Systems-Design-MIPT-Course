function main(runParam)

if runParam.runOldMusic
    perform(runParam);
end

if ~runParam.runNewMusic
    return;
end

% sample CSI trace is a 90x1 vector where first 30 elements correspond to subcarriers for first rx antenna, second 30 correspond to CSI from next 30 subcarriers and so on.
% replace sample_csi_trace with CSI from Intel 5300 converted to 90x1 vector

if runParam.isCSIfromWinner == 1
    sample_csi_traceTmp = load('vectorCSI');
    sample_csi_trace = sample_csi_traceTmp.vectorCSI;
    
    a=load('channelData.mat');
    path_power = a.P;
    aoas = -a.aoas;
    mcsDelays = a.delays * 1e6;
    arrayAngle = a.arrayAngle;
    phi = a.phi;

    if runParam.printArrayAngle
        arrayAngle
    end
    
    if runParam.printOutputMatrix
        DelayAnglePower = zeros(3, size(path_power, 2));
        
        DelayAnglePower(1,:)=mcsDelays;
        DelayAnglePower(2,:)=aoas(:,:,1);
        DelayAnglePower(3,:)=path_power;
        
        DelayAnglePower
    end
    
    if runParam.printOutput
        phi
        path_power
        mcsDelays
        aoas
    end
else
    sample_csi_traceTmp = load('sample_csi_trace');
    sample_csi_trace = sample_csi_traceTmp.sample_csi_trace;
end

if runParam.isPlotCSI == 1
   figure, plot(abs(sample_csi_trace)), title('CSI')
end

fc = runParam.fc; % center frequency
M = runParam.antennaNum;    % number of rx antennas
fs = 40e6; % channel bandwidth
c = 299792458;  % speed of light
d = c/fc/2;  % distance between adjacent antennas in the linear antenna array
% dTx = 2.6e-2;

switch runParam.subcarrierMode
    case 0
        SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,-2,2,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
    case 1
        SubCarrInd = [-58,-54,-50,-46,-42,-38,-34,-30,-26,-22,-18,-14,-10,-6,6,10,14,18,22,26,30,34,38,42,46,50,54,58]; % WiFi subcarrier indices at which CSI is available
    case 2
        SubCarrInd = [-58:-2 2:58];
    case 3
        SubCarrInd = (-58:58);
        SubCarrInd = [SubCarrInd(1:26) SubCarrInd(28:53) SubCarrInd(65:90) SubCarrInd(92:end)];
    case 4
        SubCarrInd = [-26:-1 1:26];
end

N = length(SubCarrInd); % number of subcarriers
% subCarrSize = 128;  % total number fo
fgap = 312.5e3; % frequency gap in Hz between successive subcarriers in WiFi
lambda = c/fc;  % wavelength

if runParam.doubleCSI
    T = runParam.packetNum * 2;
else
    T = runParam.packetNum;
end

% MUSIC algorithm requires estimating MUSIC spectrum in a grid. paramRange captures parameters for this grid
% For the following example, MUSIC spectrum is caluclated for 101 ToF (Time of flight) values spaced equally between -25 ns and 25 ns. MUSIC spectrum is calculated for for 101 AoA (Angle of Arrival) values between -90 and 90 degrees.
paramRange = struct;
paramRange.GridPts = runParam.grid; % number of grid points in the format [number of grid points for ToF (Time of flight), number of grid points for angle of arrival (AoA), 1]
%paramRange.delayRange = [-50 50]*1e-9; % lowest and highest values to consider for ToF grid. 
paramRange.delayRange = runParam.delayRange; % lowest and highest values to consider for ToF grid. 
paramRange.angleRange = runParam.angleRange; % lowest and values to consider for AoA grid.
do_second_iter = 0; % ???????
% paramRange.seconditerGridPts = [1 51 21 21];
paramRange.K = floor(M/2)+1; % parameter related to smoothing.
paramRange.L = floor(N/2); % parameter related to smoothing.  

if runParam.doubleCSI
    paramRange.T = runParam.packetNum * 2;
else
    paramRange.T = runParam.packetNum;
end

paramRange.deltaRange = [0 0]; 

maxRapIters = Inf;
useNoise = 0;
paramRange.generateAtot = 2;

% ToF sanitization code (Algorithm 1 in SpotFi paper)
sample_csi_trace_sanitized = zeros(size(sample_csi_trace));
for iT=1:T
    csi_plot = reshape(sample_csi_trace(:,iT), N, M);
    
    [PhsSlope, PhsCons] = removePhsSlope(csi_plot,M,SubCarrInd,N);
    ToMult = exp(1i* (-PhsSlope*repmat(SubCarrInd(:),1,M) - PhsCons*ones(N,M) ));
    csi_plot = csi_plot.*ToMult;
    relChannel_noSlope = reshape(csi_plot, N, M);
    sample_csi_trace_sanitized(:,iT) = relChannel_noSlope(:);
    
    if runParam.isPlotSanitizedCSI == 1
        figure, plot(abs(sample_csi_trace_sanitized)), title('Sanitized CSI')
    end
end

% MUSIC algorithm for estimating angle of arrival
% aoaEstimateMatrix is (nComps x 5) matrix where nComps is the number of paths in the environment. First column is ToF in ns and second column is AoA in degrees as defined in SpotFi paper
aoaEstimateMatrix = backscatterEstimationMusic(sample_csi_trace_sanitized, M, N, c, fc,...
                    T, fgap, SubCarrInd, d, paramRange, maxRapIters, useNoise, do_second_iter, runParam, ones(2)) ;
                
if runParam.printMusicOutput
    musicDelayAnglePower = aoaEstimateMatrix(:,1:3)'
end

tofEstimate = aoaEstimateMatrix(:,1); % ToF in nanoseconds
aoaEstomate = aoaEstimateMatrix(:,2); % AoA in degrees

end