clc;
clear;
close all;
yalmip('clear')

%% Load simple network
load('benchmarkNetworkSimple.mat')

%% Generate traffic matrix
cDistance = inf;
freqMax = 40;
subcarrierBandwidth = 0.0625;
totalTraffic = 200;
p_Lambda = CreateTrafficVector(networkCostMatrix, cDistance, ...
    totalTraffic, subcarrierBandwidth);
N_Lambda = size(p_Lambda, 2);

%% Run optimization - 1. preprocessing
load('benchmarkNetworkSimple.mat')
load('fiberParameter.mat')

% The bandwidth of subcarrier, unit 100 GHz
subcarrierBandwidth = 0.0625;
% The max number of subcarriers per link
subcarrierMax = round(freqMax/subcarrierBandwidth);
% The number of nodes in the network
nodeNum = size(networkAdjacentMatrix, 1);
% A vector of the spectral efficiency of the available modulation formats
% They are PM-BPSK, PM-QPSK, PM-8QAM, and PM-16QAM
sev = [2; 4; 6; 8];
demandTotalPair = 30; %nchoosek(nodeNum, 2);
% The intensity of the traffic demands
trafficDemandIntensity = [0.5, 4];
% Randomly choose demandNum pairs from the nchoosek(n£¬ 2) pairs
% vin is v_{i,n}
% 1 if n is source or destination of connection i
% size of vin: demandTotalPair x nodeNum
[demandPairs, demandPairMatrix] = createTrafficDemands(nodeNum, demandTotalPair, ...
    trafficDemandIntensity);
save('trafficMultiload_30.mat')