function [slope, lener, scales ] = waveletspectraPD(data, L, W, k1, k2, ismean, isplot)
%
%  [slope, levels, log2spec ] = waveletspectra(data, L, wf, k1, k2)
%  input:  data - data in time domain
%          L - coarse level.
%          wf - wavelet filter
%          k1 -  start with coarse level k1 when calculating slope, k1 >= L.
%          k2 -  end with the level  k2 when calculating slope, k2<=log2(n)-1
%          ismean -- 0 for median and 1 for mean, default is mean.
%  output: slope - scaling slope of log2-energies.
%          levels - integers L, L+1, ..., log2(n)-1
%          log2spec - log2 of levelwise averages of squared wavelet
%                     coefficients 
%
% 
if nargin == 1,  L=1;  wf=[sqrt(2)/2 sqrt(2)/2];  k1=1; k2=log2(length(data))-1;  ismean=1; end
if nargin == 2,        wf=[sqrt(2)/2 sqrt(2)/2];  k1=L; k2=log2(length(data))-1 ; ismean=1; end
if nargin == 3,                                   k1=L; k2=log2(length(data))-1 ; ismean=1; end
if nargin == 4,                                         k2=log2(length(data))-1 ; ismean=1; end
if nargin == 5,                                                                   ismean=1; end


n = size(data,1);
J=floor(log2(n)); 

%W=WavmatND(wf, n, J - L,  0);
wddata = W*data;
%sqttima = wdddata.^2;

lener = [];
for u= L:J-1
    ind = linspace(1,n,n)+n*(u);
    temp = wddata(ind);
    
     if  ismean == 1
        lener = [lener   log2( mean(temp.^2 )) ];
     elseif ismean == 0
        lener = [lener   log2( median(temp.^2 )) ];
        else error('not known average of energies. Use mean (ismean=1) or median (ismean=0)')
     end  
end

 slscales = k1:k2;
 sllener =  lener(round(k1-L+1):round(k2-L+1));
 b = polyfit(slscales, sllener,1);
 slope = b(1);
 intercept=b(2);
 scales = L:(J-1);
   if isplot==1
        lw = 2;
        set(0, 'DefaultAxesFontSize', 15);
        msize = 6; fs=15;
          plot(scales,lener,'bo-','linewidth',lw)
          hold on
          plot(scales,  slope*scales + intercept + 1 , 'k:','linewidth',lw)
          plot(slscales,  slope*slscales + intercept + 1, 'r-','linewidth',lw+1)
          text( k1,slope*k1+intercept+1.5, num2cell(slope) ,'fontsize',fs)
                xlabel('dyadic level','fontweight','bold','fontsize',fs)
                ylabel('log spectrum','fontweight','bold','fontsize',fs)
                axis([L-1,J, min(lener-2), max(lener+2)])
          hold off
   end 
  %else('plot not requested')
       
%-------------- Brani 10/06-------------------------------------------