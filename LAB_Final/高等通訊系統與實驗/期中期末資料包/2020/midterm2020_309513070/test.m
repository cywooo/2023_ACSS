clear;close all;clc;

x = zeros(1,32);
x(8:84)=1;
x_f = fft(x);
figure();stem(fftshift(abs(x_f)));
y = x_f((abs(x_f)<1)==0);
figure();stem(abs(y));
figure();plot(ifft(y));
