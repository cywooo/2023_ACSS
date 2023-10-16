clear all; close all; clc;
load('G_data.mat');
ph = phase(xf2);
delay = -floor((ph(end)-ph(1))/2/pi);

% rect function
x = zeros(1,511);
pos = 100;
x(pos) = 1;
% filt = ifft(xf2);
% y = conv(x,filt);
x_fft = fft(x);
y_fft = x_fft.*xf2;
y = ifft(y_fft);

figure();plot(x);hold on;plot(y);
[val,pos_d] = find(y>0.1);
xtime2 = pos_d(1)-pos;

