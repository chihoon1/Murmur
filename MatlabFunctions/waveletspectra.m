function [slope, levels, log2spec ] = waveletspectra(data, L, wf, k1, k2, ismean, isplot)
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
if nargin == 1,  L=1;  wf=[sqrt(2)/2 sqrt(2)/2];  k1=1; k2=log2(length(data))-1;  ismean=1; isplot =1; end
if nargin == 2,        wf=[sqrt(2)/2 sqrt(2)/2];  k1=L; k2=log2(length(data))-1 ; ismean=1; isplot =1; end
if nargin == 3,                                   k1=L; k2=log2(length(data))-1 ; ismean=1; isplot =1; end
if nargin == 4,                                         k2=log2(length(data))-1 ; ismean=1; isplot =1; end
if nargin == 5,                                                                   ismean=1; isplot =1; end
if nargin == 6,                                                                             isplot =1; end
lnn = floor(log2(length(data)));
%wddata = FWT_PO(data, L, wf); 
wddata = dwtr(data, lnn - L, wf); 

y = [];
for i =  L:(lnn-1)
   help = wddata(   round(2^(i)+1) : round(2^(i+1))   )   ;
   if  ismean == 1
       y = [y mean(help.^2)];
   elseif ismean == 0
       y = [y median(help.^2)];
   elseif ismean == 2
       y = [y fastDcov(help', help')];
   else error('not known average of energies. Use mean (ismean=1) or median (ismean=0)')
   end
end
   levels = L:(lnn-1);
   log2spec = log2(y);
   yy = log2spec(round(k1-L+1):round(k2-L+1));
   aa =polyfit([k1:k2], yy, 1);
   slope = aa(1);
   cc = polyval(aa,[k1:k2]);
       %--- set plotting parameters -------
       if isplot == 1
        %figure;
        lw = 2;  
        set(0, 'DefaultAxesFontSize', 15);
        fs = 15;
        msize = 6;
        %plot(levels, log2spec,  'linewidth', lw)
        hold on
        plot(levels, log2spec, 'o', 'markersize', msize)
        plot(k1:k2, yy, 'r-', 'linewidth', lw)
        plot(k1:k2, cc + 1,'r:','linewidth', lw)
        text( k1,cc(1)+0.5, num2cell(slope) )
        xlabel('dyadic level','fontweight','bold','fontsize',fs)
        %xlabel('Scales','fontweight','bold','fontsize',fs)
        ylabel('log spectrum','fontweight','bold','fontsize',fs)
        %ylabel('log energy ','fontweight','bold','fontsize',fs)
        axis tight
        hold off
       end 
       
%-------------- Brani 10/06-------------------------------------------