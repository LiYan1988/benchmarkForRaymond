function snr = standardGN(psd, f, bd, Nspan, systemParameters)
% calculate right hand side of (3) in the JLT paper
% at the input: psd is mW/THz, but in calculation psd is W/Hz

mu = systemParameters.mu;
rho = systemParameters.rho;
Nase = systemParameters.Nase;
Nuser = length(psd);

% h1
n_ase = zeros(Nuser, 1);
for i=1:Nuser
    n_ase(i) = Nase/psd(i);
end
n_ase = n_ase*1e15;

% h2
% n_sci = zeros(Nuser, 1);
% load('sciFitCoef_K30.mat', 'coefSci')
% for i=1:Nuser
%     n_sci(i) = mu*rho/pi*4*fcn_Fmm(0, 0, bd(i)*1e9, bd(i)*1e9, rho)*...
%         max([psd(i), 1]*coefSci)*1e-30;
% end

% h2 and h3
n_nli = zeros(Nuser, Nuser); % NLI nosie from j to i
for i=1:Nuser
    for j=1:Nuser
        if i==j
            n_nli(i, i) = mu*psd(i)^2*asinh(rho*(bd(i)*1e9)^2)*1e-30;
        else
            tmp = 2*abs(f(i)-f(j))/bd(j);
            n_nli(i, j) = mu*psd(j)^2*log((tmp+1)/(tmp-1))*1e-30;
        end
    end
end

n_nlitot = sum(n_nli, 2);
n_noisetot = n_ase+n_nlitot;

snr = 1./n_noisetot/Nspan;