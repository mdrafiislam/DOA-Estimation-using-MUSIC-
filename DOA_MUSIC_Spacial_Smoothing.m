%DOA estimation based on MUSIC using spatial smoothing

clc;close all;clear all;

%% Model

N=9; % # of antennas
L=5000; % # of symbols
EbN0_dB=10;
d=0.5; % Array spacing and signal wavelength ratio;
Ia1 = 30; % 1st Signal's Impringing angle
Ia2 = 45; % 2nd Signal's Impringing angle
alpha=1; % signal level of x2 relative to x1
thrs=5000; %Threshold for finding smaller eigen values
Ns=5; % # of elments in each subarray

As=N-Ns+1; % # of subarrays
RxxT=0;

%% Spatial Smoothing

for ii=0:As
    a1=exp(-1j*2*pi*([ii:ii+Ns-1].')*d*sind(Ia1));
    a2=exp(-1j*2*pi*([ii:ii+Ns-1].')*d*sind(Ia2));
    s1=((randi([0 1],1,L)*2-1) +1j*(randi([0 1],1,L)*2-1))/sqrt(2);
    x1 = a1*s1; x2 = a2*s1 * alpha;
    n = (randn(Ns,L)+1j*randn(Ns,L))/sqrt(2)* 10^(-EbN0_dB/20); % noise
    x = x1+x2+n;
    Rxxs=x*x';
    RxxT = RxxT+Rxxs;
end
Rxx = RxxT/Ns;

%% Eigen-decomposition of the covariance matrix

[U A UH]=eig(Rxx);
A_diag=diag(A);
count=1;
for ii=1:length(A_diag)
    if((A_diag(ii+1)-A_diag(ii))>thrs)
        break;
    end
    count=count+1;
end

%% Noise Subspace

Un(:,[1:count])=U(:,[1:count]); 
An(:,[1:count])=A(:,[1:count]);

angle_range=[-90:90];
for ii=1:length(angle_range)
    %a = steer(Ns,0.5,angle_range(ii));
    a=exp(-1j*2*pi*([0:Ns-1].')*d*sind(angle_range(ii)));
    p(ii) = 10*log10(abs(1/(((a'*Un)*Un')*a)));
end
figure(1)
plot(angle_range,p/max(p))
title('Array Pattern')
xlabel('Angle of arrival (degree)')
ylabel('Normalized array pattern (dB)')