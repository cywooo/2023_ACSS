clear all; clc ;close all;

n = linspace(-10,10,1001);
w = linspace(-1/2,1/2,length(n));
%% Demo 2-4
s = cos(2*pi*1*n)+3*cos(2*pi*20*n);
figure;
plot(n,s);

S = fftshift(fft(s));
figure(1);
subplot(2,1,1);
plot(w,abs(S));
title('|S(w)|');


h=zeros(1,length(n));
for i = 3:length(n)
    h(i)=1.6*h(i-1)-0.73*h(i-2)+1*(n(i)==0)+0.9*(n(i)==1);    
end

figure(2);
subplot(3,1,1);
plot(n,s);
title('s(n)');

figure(2);
subplot(3,1,2);
plot(n,h);
title('h(n)');

[conv_y,out_length] = conv_(s,h);
conv_n = linspace(2*min(n),2*max(n),out_length);
figure(2);
subplot(3,1,3);
plot(conv_n,conv_y);

S_filter = fftshift(fft(conv_y(500:1500)));
figure(1);
subplot(2,1,2);
plot(w,abs(S_filter));
title("|S'(w)|");

