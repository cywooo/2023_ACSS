clear all; close all; clc;
load('G_data.mat');

f = linspace(-1/2,1/2,length(x1));
x1_fft = fft(x1);
figure();stem(f,abs(fftshift(x1_fft)));
[u,v] = find(abs(x1_fft)>0.318);
freq1 = v/length(x1_fft);
