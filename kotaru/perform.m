function perform(runParam)

fs = 40e6;

% [N,Fo,Ao,W]=firpmord([20e3 60e3],[1 0],[(10^(0.4/20)-1)/(10^(0.4/20)+1) 10^(-200/20)],fs);
% B1=firpm(N+2,Fo,Ao,W);
% [N,Fo,Ao,W]=firpmord([40e3 80e3 150e3 200e3],[0 1 0],[10^(-200/20) (10^(0.4/20)-1)/(10^(0.4/20)+1) 10^(-200/20)], fs);
% B2=firpm(N+2,Fo,Ao,W);

B1=[1 0];
B2=[1 0];
A1=[1 0];
A2=[1 0];

L = runParam.sampleNum; % num of used samples
fc = runParam.fc;
threshold = runParam.threshold;

s=load('signal.mat');
m = s.signal';
m = m(:,2000:2000-1+L);
ind = s.fieldInd;
c = 299792458;

array = zeros(runParam.antennaNum,3);
for i=1:runParam.antennaNum
    array(i,1)=c/(2*fc)*i;
end

%[L direct dirphase dphcali m mphase mphcali]=phasesyn(m,B1,A1,B2,A2,L);   % эта штука всё портит
% из .vi от NI:
% re(dirphase) - direct wave phase before phase synchronized
% re(dphcali) - direct wave after phase synchronized
% re(mphcali) - signal phase after phase synchronized
% re(mphase) - signal phase before phase synchronized
% re(direct) - direct wave phase after phase synchronized
% re(m) - air signal after phase synchronized

[D M az el Z]=reallocation(array,fc,m,L,threshold);
% az - azimuth, el - elevation, Z - pseudospectrum, M - number of sources
% D - eigenvalues

az = az-90;

if runParam.plotOldMusic
    figure, plot(az,Z), xlabel('Angle, degrees'), ylabel('Amplitude'), title('Pseudospectrum classical Music'),
    xlim([-90 90]), ylim([min(Z)-2 max(Z)+3]), xticks(-90:30:90);
end

if runParam.printOldMusicOutput
    [Y, I] = max(Z);
    M
    peak = az(I)
end

end

%MuSIC algorithm
% x - 2D complex (the same m), L - int32
function [D M az el Z]=reallocation(array,fc,x,L,threshold)
L=double(L);
%Correlation Matrix.
Rxx=x*x'/L;
%Eigen value decomposition. E matrix with columns as eigenvectors of Rxx.
%Diagonal elements of D are eigenvalues of Rxx.
%Matrix E columns are the corresponding right eigenvectors
[E D]=eig(Rxx);
%Diagonal elements are the eigenvalues.
eigen=diag(D);
ne=[];
%Construct vectors of strong eigenvalues.
%N - antennas number (in our case 4)
N=size(array,1);

for i=1:N
    if real(eigen(i))<threshold
        ne=[ne i];
    end
end

%Number of sources (M) is equal to number of strong eigenvalues.
M=N-length(ne);
%Collect eigenvectors of the noise subspace which correspond to weak eigenvalues
En=E(:,ne);
%Vector of angle of arrivals with 0.5 deg resolution.
az=0:0.5:180;
az=az';%transpose
el=zeros(size(az));%vector of zeros with the same size as az
%Calculation of Azimuth(az) and Elevation/Altitude(el)
S=spv(array,[az,el],fc);
Z=zeros(1,length(az));
for i=1:length(az)
Z(i)=-10*log10(abs(S(:,i)'*(En*En')*S(:,i)));
end

L_out=int32(L);
end

%Make the phases of receiving boards synchronized
function [L direct2 dirphase dphcali m2 mphase mphcali]=phasesyn(m,B1,A1,B2,A2,L)
cut=1000;
direct=filter(B1,A1,m.').';%Lowpass;Get the direct wave
direct(:,L-cut+1:L)=[];
direct(:,1:cut)=[];
N=size(m,1);
dirphase=angle(direct)/pi*180;%Direct wave phase
xc=ones(N,L-2*cut);
direct2=zeros(N,L-2*cut);
direct2(1,:)=direct(1,:);
for n=2:N
    xc(n,:)=direct(1,:)./direct(n,:);
    %Phase difference between the reference receiving board and other ones
    direct2(n,:)=exp(1i*angle(xc(n,:))).*direct(n,:);
    %Direct wave after phase synchronized
end
dphcali=angle(direct2)/pi*180;%Direct wave phase after phase synchronized
m=filter(B2,A2,m.').';%Highpass;Get the signal through the air
m(:,L-cut+1:L)=[];
m(:,1:cut)=[]; 
mphase=angle(m)/pi*180;%Signal phase
m2=m;
for n=1:N 
    m2(n,:)=m2(n,:).*exp(-1i*angle(direct(n,:)));%.*abs(xc(n,:));
    %Phase synchronized
end

mphcali=angle(m2)/pi*180;%Signal phase after phase synchronised
L=L-2*cut;%The length of signal
end

function S=spv(array,direction,fc)
% Degree to radian conversion.
SOURCES=frad(direction);
% https://en.wikipedia.org/wiki/Wavenumber
% Conversion of radians to half wavelength wavenumber
% frequency * wavelength = speed of light(3e8)
%wavenumber is the number of waves in a unit distance.
%The number of radians per unit distance, sometimes termed the angular wavenumber 
%or circular wavenumber, but more often simply wavenumber (2*pi/wavelength)
%here we are calculating in terms of half wavelength.
KI = 2*fc/3e8*fki(SOURCES(:,1),SOURCES(:,2));
%array is a vector of antenna distances.
%e.g. [0 0.5 * lambda lambda 1.5 * lambda]
S = exp(-1j*(array*KI));
end

function x=frad(degrees)
%degrees to rad
x=degrees*pi/180.0;
end

function KI=fki(az,el)
%https://en.wikipedia.org/wiki/Wavenumber
%wavenumber vector in half wavelengths;
%wavenumber is the number of waves in a unit distance.
%The number of radians per unit distance, sometimes termed the angular wavenumber 
%or circular wavenumber, but more often simply wavenumber (2*pi/wavelength)
%here we are calculating in terms of half wavelength.
KI = pi*[cos(az).*cos(el)   sin(az).*cos(el)    sin(el)]';
end






