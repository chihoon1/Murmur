
close all 
clear all
clc
addpath('/Users/dixon/Documents/TAMU/DemosNew')

% [data, fs] = audioread('2530_AV.wav'); save  2530_AV.txt  data -ASCII
% [data, fs] = audioread('2530_MV.wav'); save  2530_MV.txt  data -ASCII
 [data, fs] = audioread('2530_MV.wav'); 


N=length(data);  LN= log2(N); NN=floor(LN);
y=data(1:2^NN);

figure(1)
plot(y,'-')
xlabel("Time(Seconds)"); ylabel("Aplitude")
title("An example Mitral Point(MV) Siganl")
grid on 
 
figure(2)
n1=1001; n2=2024;
yy=y(n1:n2);
plot(yy)

figure(3)
%[slope, levels, log2spec ] = waveletspectra(yy, 1, [sqrt(2)/2 sqrt(2)/2], 5,9);
subplot(221)
[slope, levels, log2spec ] = waveletspectra(y, 1, [sqrt(2)/2 sqrt(2)/2], NN-4, NN-1);
grid on 
title('MV') 

subplot(222)
[data, fs] = audioread('2530_AV.wav'); 
N=length(data);  LN= log2(N); NN=floor(LN);
y=data(1:2^NN);

[slope, levels, log2spec ] = waveletspectra(y, 1, [sqrt(2)/2 sqrt(2)/2], NN-4, NN-1);
grid on
title('AV')

subplot(223)
[data, fs] = audioread('2530_PV.wav'); 
N=length(data);  LN= log2(N); NN=floor(LN);
y=data(1:2^NN);

[slope, levels, log2spec ] = waveletspectra(y, 1, [sqrt(2)/2 sqrt(2)/2], NN-4, NN-1);
grid on
title('PV')

subplot(224)
[data, fs] = audioread('2530_TV.wav'); 
N=length(data);  LN= log2(N); NN=floor(LN);
y=data(1:2^NN);

[slope, levels, log2spec ] = waveletspectra(y, 1, [sqrt(2)/2 sqrt(2)/2], NN-4, NN-1);
grid on
title('TV')


