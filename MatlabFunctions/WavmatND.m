function W = WavmatND(hf, m, p, shift)
% WavmatND -- Transformation Matrix for Non-Decimated WT
%  Usage
%    dat=[1 0 -3 2 1 0 1 2];
%    wavfilt = [1/sqrt(2) 1/sqrt(2)];
%    W = WavmatND(wavfilt,8,3,2);
%    wt = W * dat'  % to obtain a transformed signal
%    data = W' * wt % should return you to the 'dat'
% 
%  Inputs
%    filter   --   wavelet filter
%    N        --   size of matrix/length of data.        
%    p        --   depth of transformation. Range >= 1. Although  not
%                  limited p > log2(N) leads to oversmoothing and
%                  high influence of the filter
%    shift    --   the matrix is not unique and any integer shift gives
%                   a valid transformation. Default is 0.
% 
%  Output
%    W        --  (p + 1)*N x N transformation matrix 
%
%  Description
%    For a quadrature mirror filter h (low pass) the wavelet
%    matrix is formed to perform Non-decimated wavelet transform.  
%
%  Copyright (C) -- see WavmatND/Copyright

    hf=hf(:)';  gf = fliplr(conj(hf).* (-1).^(0:length(hf)-1 ));
    W=[]; hmatold=eye(m);
    h=[hf,zeros(1,m)]; %extend filter H by 0's to sample by modulus
    g=[gf,zeros(1,m)]; %extend filter G by 0's to sample by modulus
for i = 1:p
    clear gmat; clear hmat;
    for  jj= 1:m
       for ii=1:m
           modulus = mod(m + ii - jj - shift , m) + 1;
           modulus = modulus + (modulus == 0)*m;   
           hmat(ii,jj) =  h(modulus);
           gmat(ii,jj) = g(modulus);
       end
    end
%-------------------------------------
    W=[ gmat' * hmatold  ; W];
    smooth = hmat'* hmatold ;
    hmatold = smooth; 
    h = [dilate_filter(hf,2^(i)-1),zeros(1,m)];
    g = [dilate_filter(gf,2^(i)-1),zeros(1,m)];
end

W=[smooth; W];
end
%--------------------------------------
function  filtd = dilate_filter(filt,k)
%-------------------------------------
%  Dilate Filter by k zeros between the original 
% taps, k integer > 0.
%-------------------------------------
newlength = (k+1)*length(filt)-k;
filtd = zeros(1,newlength);
filtd(1:(k+1):newlength)=filt;
%--------------------------------------
% 
end
function out=imbed(a)
% imbeds 0 in the filter a
b=[a; repeat(0, length(a))];
out=b(:)';
end
function b = repeat(a, n)
% Repeats an array a n times
% Usage
%   b = repeat(a, n)
% Input
%   a, n
% Output
%   b
%  
b=[];
for i=1:n
b=[b a];
end
end
