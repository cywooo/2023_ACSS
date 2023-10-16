clear all
close all

BT = 0.5;
M = 16;
trun = 5;
h = Gaussian_filter(BT,M,trun);
figure(1);
plot(h);
title('Impulse response');

figure(2);
plot(fftshift(abs(fft(h))));
title('Frequency response');