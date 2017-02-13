% calculate GN model with linearized model
clc;
clear;
close all;

%% define fiber parameters
alpha = 0.22; % dB/km, attenuation of fiber, NOTE: alpha is positive!
alpha = alpha*1e-4*log(10); % 1/m
L = 100e3; % m, length of one span
h = 6.626e-34; % J*s, Plank's constant
niu = 193.548e12; % Hz, frequency of lightwave at 1550 nm
nsp = 10^(5.5/10)/2; % spontaneous emission factor
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

%% test linearized GN model
Nuser = 10;
psd = randi([10, 20], [Nuser, 1]);
f = (0:1:Nuser-1)*200;
bd = 50*ones(Nuser, 1);
Nspan = 10;
snrL = psd./linearizedGN(psd, f, bd, Nspan, systemParameters);

%% test standard GN model
snrS = psd./standardGN(psd, f, bd, Nspan, systemParameters);

figure; hold on;
plot(snrL)
plot(snrS, '--')