clear all; clc; close all;
%% practice 1
n = [1:201];
w = linspace(-1/2,1/2,length(n));
x = cos(2*pi*1/64*n);
carrier = cos(2*pi*1/4*n);

y = x.*carrier;

figure(1);
subplot(6,1,1);
plot(n,x);
figure(1);
subplot(6,1,2);
plot(n,carrier);
figure(1);
subplot(6,1,3);
plot(n,y);

y_receive = y.*carrier;
X = fftshift(fft(x));
Y = fftshift(fft(y));
Y_receive = fftshift(fft(y_receive));

figure(2);
subplot(5,1,1);
plot(w,abs(X));
figure(2);
subplot(5,1,2);
plot(w,abs(Y));
figure(2);
subplot(5,1,3);
plot(w,abs(Y_receive));

LPF = 0 + 2 .*(w > -0.2).*(w < 0.2);
Y_de = Y_receive.*LPF;

figure(2);
subplot(5,1,4);
plot(w,abs(LPF));
figure(2);
subplot(5,1,5);
plot(w,abs(Y_de));

y_de = ifft(fftshift(Y_de));
y_receive = ifft(fftshift(Y_receive));

figure(1);
subplot(6,1,4);
plot(n,abs(x));
figure(1);
subplot(6,1,5);
plot(n,abs(y_receive));
figure(1);
subplot(6,1,6);
plot(n,abs(y_de));

