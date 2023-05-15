function H = wavelet_entropy(x, filt, L)
    J = log2(size(x,1));
    % Compute the wavelet transform of the signal
    wt_x = dwtr(x, J - L, filt);
    %C = wt_x(2^(J-L)+1: 2^((J-L)+1));
    E = [];
    for i =  L:(J-1)
        C = wt_x( round(2^(i)+1) : round(2^(i+1))   );
        % Compute the energy of the wavelet coefficients
        E = [E C.^2];
        % Normalize the energy
    end 
    P = E/sum(E);
    P(P==0) = [];
    % Compute the entropy
    H = -sum(P.*log2(P));
end