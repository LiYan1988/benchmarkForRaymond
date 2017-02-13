function result = fcn_bm(freqMax, demandPairs, demandPairMatrix, psdfix)
% variable PSD
yalmip('clear')

%% Load network data
load('GermenNetworkTopology.mat')
load('fiberParameter.mat')

% The max number of subcarriers per link
subcarrierMax = freqMax;

% The number of nodes in the network
% nodeNum = size(networkAdjacentMatrix, 1);

% A vector of the spectral efficiency of the available modulation formats
% They are BPSK, PM-QPSK, PM-8QAM, and PM-16QAM
sev = [1];
% number of total demands
demandTotalPair = size(demandPairs, 1); 

% The intensity of the traffic demands
% trafficDemandIntensity = [0.5, 4];

% Randomly choose demandNum pairs from the nchoosek(n£¬ 2) pairs
% vin is v_{i,n}
% 1 if n is source or destination of connection i
% size of vin: demandTotalPair x nodeNum
% [~, demandPairMatrix] = createTrafficDemands(nodeNum, demandTotalPair, ...
%     trafficDemandIntensity);
p_Lambda = demandPairs(:, 3);

%% Create z_{l,n}, p_incidence
% 1 if node n is an ending node of link 
% size: linkNum x nodeNum
% zln is the incidence matrix of an bidirectional network, all entries are
% 1.
[~, edgeIdxList, p_incidence, linkLength] = ...
    findEdge(networkAdjacentMatrix, networkCostMatrix);
linkLength = linkLength/100;
linkNum = length(edgeIdxList);

%% Create variable pil
% 1 if link l is on the route of connection i
% size: demandTotalPair x linkNum
x_connectionLink = binvar(demandTotalPair, linkNum);

%% Create variable qin
% 1 if node n is on the route of connection i
% size: demandTotalPair x nodeNum
% x_connectionNode = binvar(demandTotalPair, nodeNum);

%% Create constraint (3e)
F = [x_connectionLink*p_incidence==demandPairMatrix];
F = [F; sum(x_connectionLink, 2)<=linkNum/2];

%% Create variable yij
% 1 if connection i and j share a link
% the index corresponding to (i, j) is
% k = (2*demandTotalPair-i)*(i-1)/2+j-i, j>i
% the actual variable is the upper triangular matrix
x_connectionShare = binvar(nchoosek(demandTotalPair, 2), 1);

%% Create constraint (3f)
for i = 1:demandTotalPair-1
    for j = i+1:demandTotalPair
        for l = 1:linkNum
            k = (2*demandTotalPair-i)*(i-1)/2+j-i;
            F = [F; x_connectionShare(k)+1>=...
                x_connectionLink(i, l)+x_connectionLink(j, l)];
        end
    end
end

%% Create variable u_{ij}
% 1 if yij=1 and fi+Bi<fj
% the index corresponding to (i, j) is
% k = (2*demandTotalPair-i)*(i-1)/2+j-i, j>i
% the actual variable is the upper triangular matrix
x_connectionOdering = binvar(nchoosek(demandTotalPair, 2), 1);

% The range of x_connectionOdering, uij is possible to be 1 only if yij is
% 1, and uij can be 0 if yij is 1.
F = [F; x_connectionOdering<=x_connectionShare];

%% Create variables fi
% the lowest frequency allocated to connection i
x_lowestFrequency = sdpvar(demandTotalPair, 1);

% The range of fi
F = [F; x_lowestFrequency>=0; ...
    x_lowestFrequency<=subcarrierMax];

%% Create variables binary for modulation format, and (3b)
% 1 if connection i is assigned modulation format k
x_modulationBinary = binvar(demandTotalPair, length(sev));

% Asscociated constraint with x_modulationBinary
for i = 1 : demandTotalPair
    F = [F; sum(x_modulationBinary(i, :))==1];
end

%% Create variable Bi, (3c)
% Bandwidth of connection i, in number of subcarriers, implicitly integer.
x_connectionBandwidth = sdpvar(demandTotalPair, 1);

% The range of x_connectionBandwidth
F = [F; x_connectionBandwidth>=p_Lambda/sev(end); ...
    x_connectionBandwidth<=p_Lambda/sev(1)];

% The choise of x_connectionBandwidth is related with the modulation
% formats.
for i = 1 : demandTotalPair
    F = [F; x_connectionBandwidth(i)==...
        x_modulationBinary(i, :)*p_Lambda(i)./sev];
end

%% Create variable x_channelSpacing
% Channel spacing in unit of number of subcarriers, 1/2-integer
x_channelSpacing = sdpvar(nchoosek(demandTotalPair, 2), 1);
ttemp = sort(p_Lambda);
% If i and j are on the same link, their lowest possible channel spacing
% is (p_Lambda(i)+p_Lambda(j))/sev(end)/2.
% But if i and j are not on the same link, their channel spacing is 0. This
% should be consistant with the calculation of XCIs.
F = [F; x_channelSpacing>=0; ...
    x_channelSpacing<=subcarrierMax-...
    sum(ttemp(1)/sev(end)+ttemp(2)/sev(end))/2+1];
clear ttemp

%% Create constraint: ordering and channel spacings of connections
% For connections i and j, they share one link iff. yij==1. If yij==1 &
% uij==1, then j is at the right side of i; if yij == 1 & uij==0, then j is
% at the left side of i.
for i = 1:demandTotalPair-1
    for j = i+1:demandTotalPair
        k = (2*demandTotalPair-i)*(i-1)/2+j-i;
        F = [F; ...
            % yij+uij==2, fi+Bi<=fj
            implies(x_connectionShare(k)+x_connectionOdering(k)==2,...
            x_lowestFrequency(i)+x_connectionBandwidth(i)<=...
            x_lowestFrequency(j));...
            % channel spacing
            implies(x_connectionShare(k)+x_connectionOdering(k)==2,...
            x_channelSpacing(k)==...
            x_lowestFrequency(j)+x_connectionBandwidth(j)/2-...
            (x_lowestFrequency(i)+x_connectionBandwidth(i)/2));
            % yij+uij==1, fj+Bj<=fi
            implies(x_connectionShare(k)+x_connectionOdering(k)==1,...
            x_lowestFrequency(j)+x_connectionBandwidth(j)<=...
            x_lowestFrequency(i));...
            % channel spacing
            implies(x_connectionShare(k)+x_connectionOdering(k)==1,...
            x_channelSpacing(k)==...
            x_lowestFrequency(i)+x_connectionBandwidth(i)/2-...
            (x_lowestFrequency(j)+x_connectionBandwidth(j)/2))];
    end
end

%% Create variables for sharing link
% x_connectionShareLink(k, l) = 1 if connection i and j share link l, and
% k = (2*demandTotalPair-i)*(i-1)/2+j-i.
x_connectionShareLink = binvar(nchoosek(demandTotalPair, 2), linkNum);
for l = 1 : linkNum
    for i = 1 : demandTotalPair-1
        for j = i+1 : demandTotalPair
            idxtemp = (2*demandTotalPair-i)*(i-1)/2+j-i;
            F = [F; x_connectionShareLink(idxtemp, l)<=x_connectionLink(i, l);...
                x_connectionShareLink(idxtemp, l)<=x_connectionLink(j, l);...
                x_connectionShareLink(idxtemp, l)>=...
                x_connectionLink(i, l)+x_connectionLink(j, l)-1];
        end
    end
end

% x_connectionShareLink(k, l) <= x_connectionShare(k)
for k = 1 : nchoosek(demandTotalPair, 2)
    F = [F; x_connectionShareLink(k, :) <= x_connectionShare(k)];
end

%% Create variable, PSD
x_connectionPSD = sdpvar(demandTotalPair, 1);
F = [F; x_connectionPSD==psdfix];

%% Create variables for XCI term
% The NLI terms, the off-diagonal entry (i, j) represents the XCI from
% connection j to i. The diagonal entry (i, i) represents the SCI of
% connection i. This matrix is not symmetric. And this is calculated on
% each link.
x_nliLink = sdpvar(demandTotalPair, demandTotalPair, linkNum, 'full');

% The range of x_nliLink
F = [F; x_nliLink(:)>=0; x_nliLink<=max(linkLength)*...
    miu*psdfix^2*max(asinh(rho*(max(p_Lambda)/sev(1))^2), log(2.05/0.05))];

% Load fitting coefficients
load('xciFitCoef_K60.mat', 'coefXci')
load('sciFitCoef_K30.mat', 'coefSci')
% load('logFitCoef.mat')
for l = 1:linkNum
    for i = 1:demandTotalPair
        for j = 1:demandTotalPair
            if i==j
                % SCI term
                % For each modulation format of i, there is an 
                % approximation for SCI
                for k = 1:length(sev)
                    bw = ceil(p_Lambda(j)/sev(k));
                    F = [F; implies(x_modulationBinary(j, k)+...
                        x_connectionLink(j, l)==2, ...
                        x_nliLink(i, j, l)>=linkLength(l)*miu*rho/pi*4*...
                        fcn_Fmm(0, 0, bw, bw, rho)*...
                        max([x_connectionPSD(j), 1]*coefSci))];
                end
            else
                % XCI term
                % For each modulation format of j, there is an 
                % approximation for XCI
                for k = 1:length(sev)
                    % Calculate the index for connection Share
                    if j>i
                        idxtemp = (2*demandTotalPair-i)*(i-1)/2+j-i;
                    else
                        idxtemp = (2*demandTotalPair-j)*(j-1)/2+i-j;
                    end
                    F = [F; implies(x_modulationBinary(j, k)+...
                        x_connectionShareLink(idxtemp, l)==2, ...
                        x_nliLink(i, j, l)>= linkLength(l)*miu*...
                        max([x_connectionPSD(j), ...
                        2*x_channelSpacing(idxtemp)/...
                        p_Lambda(j)/sev(k), 1]*coefXci))];
                end
            end
        end
    end
end

% Create variables for NLI of all links
x_nli = sdpvar(demandTotalPair, demandTotalPair, 'full');
for i = 1 : demandTotalPair
    for j = 1 : demandTotalPair
        F = [F; x_nli(i, j)==sum(x_nliLink(i, j, :))];
    end
end

%% Create variables for ASE term
load('aseFitCoef_K30.mat', 'coefAse')
x_aseLink = sdpvar(demandTotalPair, linkNum);
F = [F; x_aseLink>=0; x_aseLink<=max(linkLength)*Nase/psdfix];
for i = 1 : demandTotalPair
    for l = 1 : linkNum
        F = [F; implies(x_connectionLink(i,l)==1,...
            x_aseLink(i, l)>=linkLength(l)*Nase*...
            max([x_connectionPSD(i), 1]*coefAse))];
    end
end

x_ase = sdpvar(demandTotalPair, 1);
for i = 1 : demandTotalPair
    F = [F; x_ase(i)==sum(x_aseLink(i, :))];
end

%% Create variables for SNR threshold term, -1/SNRth
load('snrthCoef.mat')
x_snrth = sdpvar(demandTotalPair, 1);
for i = 1 : demandTotalPair
    F = [F; x_snrth(i)==x_modulationBinary(i, :)*(snrthCoef(1:length(sev)).')];
end

%% Create variables for the whole SNR constraint
x_snrcon = sdpvar(demandTotalPair, 1);
for i = 1 : demandTotalPair
    F = [F; x_snrcon(i) == sum(x_nli(i, :))+x_ase(i)+x_snrth(i)];
end
F = [F; x_snrcon <= 0];

%% Create transmission reach constraint
% If mik==1, i.e., connection i uses modulation k, the number of total
% spans should be less than a certain number.
for i = 1 : demandTotalPair
    for k = 1 : length(sev)
        tempNoise = 3*(miu*asinh(rho*(p_Lambda(i)/sev(k))^2)*Nase^2/4)^(1/3);
        tempSpanLimit = floor(-snrthCoef(k)/tempNoise);
        F = [F; implies(x_modulationBinary(i, k)==1, ...
            x_connectionLink(i, :)*linkLength<=tempSpanLimit)];
    end
end

%% Create variable x_totalBandwidth
x_totalBandwidth = sdpvar(1);
F = [F; x_totalBandwidth<=subcarrierMax];

%% Estimation of lower bound
% The total bandwidth should be greater than the sum of all the channels'
% bandwidth link
x_connectionModulationLink = binvar(demandTotalPair, length(sev), linkNum);
for i = 1 : demandTotalPair
    for k = 1 : length(sev)
        for l = 1 : linkNum
            F = [F; x_connectionModulationLink(i, k, l)<=...
                x_modulationBinary(i, k);...
                x_connectionModulationLink(i, k, l)<=...
                x_connectionLink(i, l);...
                x_connectionModulationLink(i, k, l)>=...
                x_modulationBinary(i, k)+x_connectionLink(i, l)-1];
        end
    end
end

bik = ceil(bsxfun(@rdivide, p_Lambda, sev.'));
x_linkTotalBandwidth = sdpvar(linkNum, 1);
for l = 1 : linkNum
    F = [F; x_linkTotalBandwidth(l)==...
        sum(sum((bik.*x_connectionModulationLink(:, :, l))))];
end
F = [F; x_totalBandwidth >= max(x_linkTotalBandwidth)];

% Define total bandwidth
F = [F; x_totalBandwidth>=x_connectionBandwidth+x_lowestFrequency-1];

%% Create
w = 1e-2; w2 = 0e-4;
Objective = x_totalBandwidth+w*(sum(x_nliLink(:))+sum(x_aseLink(:)))+...
    w2*sum(x_linkTotalBandwidth);

%% Run yalmip
options = sdpsettings('solver', 'gurobi', 'gurobi.symmetry', 1, ...
    'gurobi.mipfocus', 1, 'gurobi.timelimit', 300, 'gurobi.MIPGapAbs', 1);
optimize(F, Objective, options);

result = struct();
result.ase = value(x_ase);
result.aseLink = value(x_aseLink);
result.channelSpacing = value(x_channelSpacing)*subcarrierBandwidth;
result.connectionBandwidth = value(x_connectionBandwidth)*subcarrierBandwidth;
result.connectionLink = value(x_connectionLink);
result.connectionModulationLink = value(x_connectionModulationLink);
% result.connectionNode = value(x_connectionNode);
result.connectionOdering = value(x_connectionOdering);
result.connectionPSD = value(x_connectionPSD);
result.connectionShare = value(x_connectionShare);
result.connectionShareLink = value(x_connectionShareLink);
result.lowestFrequency = (value(x_lowestFrequency)-1)*subcarrierBandwidth;
result.modulationBinary = value(x_modulationBinary);
result.nli = value(x_nli);
result.nliLink = value(x_nliLink);
result.snrcon = value(x_snrcon);
result.snrth = value(x_snrth);
result.totalBandwidth = value(x_totalBandwidth)*subcarrierBandwidth;
end