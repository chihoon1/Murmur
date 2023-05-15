function [scales,  lener,  slope] = WaveletSpectra2ND(x, L, qmf, k1, k2, isreal, isplot)
% WaveletSpectra2NDM -- 2-D  non-decimated scale-mixing wavelet spectra
%  Usage
%    [scales  logene] = WaveletSpectra2NDM(x,L,qmf,k1,k2, trim, isplot)
% 
%  Inputs
%    x   -- 2-d image (n by n real array, n dyadic = 2^J)
%    L   -- degree of coarsest scale
%    qmf -- orthonormal quadrature mirror filter 
%    k1  -- minimal level for slope calculation (default L)
%    k2  -- maximal level for slope calculation (default J-1)
%    trim -- trimmed mean parameter (default 0, plain mean)
%    isreal -- imaginary/real/ phase/magnitude of the complex wt
%    isplot -- 1 if plot;  0 if not (default 1)
% 
%  Outputs
%    scales -- scales (L,J-1) 
%    lener  -- logarithm of average energies
%    slope  -- regression fit results
% 
% Copyright (C) -- see WavmatND/Copyright

[m, n]=size(x);   J=floor(log2(min(m,n))); 

if nargin==5;  isreal = 1;isplot=1; end
if nargin==4;  k2 = J-1; isreal = 1;  isplot=1; end
if nargin==3;  k1 = L;   k2 = J-1;    isreal = 1; isplot=1; end
if nargin==2;  qmf = [sqrt(2)/2 sqrt(2)/2]; k1=L; k2=J-1;isreal = 1;  isplot=1; end
if nargin==1;  L=1; qmf=[sqrt(2)/2 sqrt(2)/2]; k1=L; k2=J-1; isreal = 1; isplot=1; end

W1=WavmatND(qmf,m,J-L,0);
W2=WavmatND(qmf,n,J-L,0);
ttima = W1*x*W2';

if isreal == 0
    wddata = imag(ttima);
elseif isreal == 1
    wddata = real(ttima);
elseif isreal == 2
    wddata  = atan(imag(ttima)./real(ttima));
elseif isreal == 3
    wddata = sqrt(imag(ttima).^2 + real(ttima).^2);
elseif isreal == 4
    wddata = sqrt(imag(ttima).^2 + real(ttima).^2) + atan(imag(ttima)./real(ttima));
elseif isreal == 5
    wddata = imag(ttima) + atan(imag(ttima)./real(ttima));
elseif isreal == 6
    wddata = real(ttima) + atan(imag(ttima)./real(ttima));

else
    error('not known real or imaginary part. Use isreal=0,1, 2, 3, 4, 5 or 6');

end 

sqttima= wddata.^2;
scales = L:(J-1);

lener=[];
for u=1:J-L
    row = linspace(1,m,m)+m*(u);
    col = linspace(1,n,n)+n*(u); 
    temp= sqttima(row,col);
    lener =[lener   log2( mean(temp(:))) ];
end
  slscales = k1:k2;
  sllener =  lener(k1:k2);
  b = polyfit(slscales, sllener,1);
  slope = b(1);
  intercept=b(2);
   if isplot==1
        lw = 2;
        set(0, 'DefaultAxesFontSize', 15);
        msize = 6; fs=15;
          plot(scales,lener,'bo-','linewidth',lw)
          hold on
          plot(scales,  slope*scales + intercept + 1 , 'k:','linewidth',lw)
          plot(slscales,  slope*slscales + intercept + 1, 'r-','linewidth',lw+1)
          text( k1,slope*k1+intercept+1.5, num2cell(slope) )
                xlabel('dyadic level','fontweight','bold','fontsize',fs)
                ylabel('log spectrum','fontweight','bold','fontsize',fs)
                axis([L-1,J, min(lener-2), max(lener+2)])
          hold off
  %else('plot not requested')
      
    end
