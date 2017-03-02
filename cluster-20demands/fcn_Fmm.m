function fmm = fcn_Fmm(cf1, cf2, bw1, bw2, rho)
    % Calculate the exact dilog function in Pontus' GN model.
    %
    % INPUT:
    % cf1, cf2: the central frequencies of the two channels
    % bw1, bw2: the bandwidth of the second channel, which generate XCI to 
    %           the first channel.
    % rho: fiber parameter, rho = xi/8, where xi is in Pontus' paper.
    %
    % OUTPUT:
    % fmm: the value of the F^{2}_{mm'} function.
    
    xi = rho*8;
    x1 = bw1.*(cf1-cf2+bw2/2)*xi/2;
    x2 = bw1.*(cf2-cf1+bw2/2)*xi/2;
    
    fmm = abs(2*(imag(dilog(1i*x1))+imag(dilog(1i*x2)))/xi);
end