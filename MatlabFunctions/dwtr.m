function dwtr = dwtr(data, L, filterh)
%  function dwtr = dwtr(data, L, filterh); 
%  Calculates the DWT of periodic data set
%  with scaling filter  filterh  and  L  detail levels. 
%
%   Example of Use:
%   data = [1 0 -3 2 1 0 1 2]; filterh = [sqrt(2)/2 sqrt(2)/2];
%   wt = dwtr(data, 3, filterh)
% 
%     1.4142   -1.4142    1.0000  -1.0000         0.7071   -3.5355
%     0.7071                  -0.7071
% 
%--------------------------------------------------------------------------

n = length(filterh);                  %Length of wavelet filter
C = data(:)';                             %Data (row vector) live in V_j
dwtr = [];                                  %At the beginning dwtr empty                
H = filterh;
G = fliplr(filterh);
G(1:2:n) = -G(1:2:n);              %    counterpart
for j = 1:L                                 %Start cascade 
   nn = length(C);                    %Length needed to 
   C = [C(mod((-(n-1):-1),nn)+1)  C  C]; % make periodic
      %D=filter(G,[1],C);                  %Filter with high-pass
   D = conv(C,G);                     %Convolve is equivalent,
   D = D([n:2:(n+nn-2)]+1);     %    keep periodic and decimate
      %C=filter(H,[1],C);                    %Filter with low-pass
   C = conv(C,H);                     %or Convolve 
   C = C([n:2:(n+nn-2)]+1);     %    keep periodic and decimate
   dwtr = [D,dwtr];                   %Add detail level to dwtr                          
end;                                          %Back to cascade or end
dwtr = [C, dwtr];                     %Add the last ``smooth'' part

