function H = wavelet_Lentropy(x, filt, L)
    J = log2(size(x,1));
    % Compute the wavelet transform of the signal
    wt_x = dwtr(x, J - L, filt);
    
    PE = [];
    for i = 1:J-L
        C = wt_x(2^(i)+1: 2^(i+1));
        % Compute the energy of the wavelet coefficients
        E = C.^2;
        % Normalize the energy
        P = E/sum(E);
        PE = [PE P];
    end 
    PE(PE==0) = [];
    %PE = PE/sum(PE);
    % Compute the entropy
    H = -sum(PE.*log2(PE));
    %H = info_entropy (PE(1:end), 'bit');
    
end