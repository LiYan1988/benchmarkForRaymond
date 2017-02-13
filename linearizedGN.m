function snr = linearizedGN(psd, f, bd, Nspan, systemParameters)
% calculate right hand side of (3) in the JLT paper
% at the input: psd is mW/THz, but in calculation psd is W/Hz

mu = systemParameters.mu;
rho = systemParameters.rho;
Nase = systemParameters.Nase;
Nuser = length(psd);

% h1
n_ase = zeros(Nuser, 1);
load('aseFitCoef_K30.mat', 'coefAse')
for i=1:Nuser
    n_ase(i) = max([psd(i), 1]*coefAse);
end
n_ase = Nase*n_ase*1e15;

% h2
% n_sci = zeros(Nuser, 1);
% load('sciFitCoef_K30.mat', 'coefSci')
% for i=1:Nuser
%     n_sci(i) = mu*rho/pi*4*fcn_Fmm(0, 0, bd(i)*1e9, bd(i)*1e9, rho)*...
%         max([psd(i), 1]*coefSci)*1e-30;
% end

% h2 and h3
n_nli = zeros(Nuser, Nuser); % NLI nosie from j to i
load('sciFitCoef_K30.mat', 'coefSci')
load('xciFitCoef_K60.mat', 'coefXci')
for i=1:Nuser
    for j=1:Nuser
        if i==j
            n_nli(i, i) = mu*rho/pi*4*fcn_Fmm(0, 0, bd(i)*1e9, bd(i)*1e9, rho)*...
                max([psd(i), 1]*coefSci)*1e-30;
        else
            n_nli(i, j) = mu*max([psd(j), 2*abs(f(i)-f(j))/bd(j), 1]*coefXci)*1e-30;
        end
    end
end

n_nlitot = sum(n_nli, 2);
n_noisetot = n_ase+n_nlitot;

snr = 1./n_noisetot/Nspan;