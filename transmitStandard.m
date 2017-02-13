% transmit 11 channels with variable bandwidth (30~350 GHz), 10 GHz
% guardband on 1 span and calculate the noise at the center channel
clc;
clear;
close all;

%% define fiber parameters
alpha = 0.22; % dB/km, attenuation of fiber, NOTE: alpha is positive!
alpha = alpha*1e-4*log(10); % 1/m
L = 100e3; % m, length of one span
h = 6.626e-34; % J*s, Plank's constant
niu = 193.548e12; % Hz, frequency of lightwave at 1550 nm
nsp = 10^(5/10)/2; % spontaneous emission factor
Nase = (exp(alpha*L)-1)*h*niu*nsp; % ASE per polarization per span
% W/Hz, signal side ASE noise spectral density
gamma = 1.32e-3; % 1/(W*m), nonlinear parameter
% gamma = 0;
beta = -2.1668e-26; % s^2/m, GVD parameter, D = 18 ps/(nm*km),
% beta = -D*lambda^2/(2*pi*c)
beta = abs(beta); % the absolute value is used in calculation

mu = 3*gamma^2/(2*pi*alpha*abs(beta));
rho = pi^2*abs(beta)/(2*alpha);

systemParameters = struct();
systemParameters.mu = mu;
systemParameters.rho = rho;
systemParameters.Nase = Nase;

%% test standard GN model
bdv = 30:10:350;
noiseCenterStandard = zeros(size(bdv));
for i=1:length(bdv)
    Nuser = 11;
    bd0 = bdv(i);
    psd = 15*ones(Nuser, 1); % mW/THz
    bd = bd0*ones(Nuser, 1);
    f = (-5:5)*(bd0+10);
    Nspan = 1;
    noiseStandard = psd./standardGN(psd, f, bd, Nspan, systemParameters)*1e-15;
    noiseCenterStandard(i) = noiseStandard(6);
    noiseLinearized = psd./linearizedGN(psd, f, bd, Nspan, systemParameters)*1e-15;
    noiseCenterLinearized(i) = noiseLinearized(6);
end

%% error
err = (noiseCenterLinearized./noiseCenterStandard-1);
figure; hold on;
plot(bdv, err)

%% 
figure;
plot(bdv, noiseCenterStandard)
