function [f,alfa]=mfspectra(data, L, wf, q, k1, k2)
%function [f,alfa]=mfspectra(data, L, wf, q, k1, k2)
%
%   Input:data---the input signal vector;
%         q------the range of moment orders (usually -1:0.1:6, or similar)
%         k1,k2-----the range of interest scale(min and max)
%         L,wf----parameters for wavlet decomposition DWTR, depth and
%         filter
%   Output: f----multifractal spectrum f(alpha)
%   For example: [a,b]=mfspectra(m,-1:0.2:6,3,10);
%             where m is the fractal signal, such as fbm(1/3) with 2^16 
%             length. q=[-1,6] with equal space 0.2. [j1, j2]=[3,10]
%   reference:P. Goncalves, R. H. Riedi and R. G. Baraniuk 
%            Simple Statistical Analysis of Wavelet-based Multifractal 
%            Spectrum Estimation, Proceedings of the 32nd Conference on 
%           `Signals, Systems and Computers', Asilomar, Nov 1998 
%--------------------------------------------------------------------------------------
%if nargin == 4,  L=5; wf = [1/sqrt(2) 1/sqrt(2)]; end;  %default
%paramteres
lnn = log2(length(data));
%wddata = FWT_PO(data, L, wf); %or
wddata = dwtr(data, lnn-L, wf); 
for i =  L:(lnn-1)
    j=i;
    help = 2^(j/2)*wddata((2^(i)+1):(2^(i+1))); %L1 normalization
    for k=1:length(q)
        s(j,k)=mean(abs(help).^q(k)); %Partition Function
    end;
end;
t=[];
for k=1:length(q)
a=polyfit(k1:k2,log2(s(k2:-1:k1,k))',1); %regression 
t=[t,a(1)];
end;
alfa=diff(t)./diff(q); % numerical derivative
qq=q(1:length(q)-1);
f=qq.*alfa-t(1:length(q)-1);
