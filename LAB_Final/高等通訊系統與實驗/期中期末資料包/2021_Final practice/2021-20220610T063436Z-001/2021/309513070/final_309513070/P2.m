clear all; close all; clc;
% load signal
load('a_2.mat');

% upsample pass SRRC filter
uf1 = 16;
x1 = zeros(1,length(a_2)*uf1);
x1(1:uf1:end) = a_2;
srrc = rcosdesign(0.25,6,16);
x2 = conv(x1,srrc);

% DMA filter
uf2 = 8;
x3 = zeros(1,length(x2)*uf2);
x3(1:uf2:end) = x2;
x4 = filter(p2filter,x3);
x4f = fftshift(fft(x4));
% figure(1);
% x3f = fftshift(fft(x3));
% f = linspace(-0.5,0.5,length(x3f));
% stem(f,abs(x3f));
% figure(2);
% f = linspace(-0.5,0.5,length(x4f));
% stem(f,abs(x4f));

% up conversion
wc = 32/128 * 2*pi;
x_uc = x4.*exp(1j*wc*[0:length(x4)-1]);
x_uc = real(x_uc);

% channel
y = x_uc;

% IF demod
wif = 2/128 * 2*pi;
y1 = y.*cos((wc-wif)*[0:length(y)-1]);
% y1f = fftshift(fft(y1));
% figure(3);
% f = linspace(-0.5,0.5,length(y1f));
% stem(f,abs(y1f));


% DMA filter
delay = 9;
y2 = filter(p2filter,y1);
% figure(4);
% y2f = fftshift(fft(y2));
% f = linspace(-0.5,0.5,length(y2f));
% stem(f,abs(y2f));
y2 = y2(1+delay:uf2:end);

% down conversion
y3 = y2.*exp(-1*1j*wif*8*[0:length(y2)-1]);

% SRRC filter
delay_sc = 96;
y = conv(y3,srrc);
a_2h = y(1+delay_sc:uf1:end);

figure();
stem(a_2);
hold on;
stem(30*a_2h(1:end));


