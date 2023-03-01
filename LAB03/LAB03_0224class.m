clear all; clc ;close all;

n =  linspace(-10,10,1001);

%% LAB 03 Practice 1 Q1
s = cos(2*pi*0.2*n);
%s = t.*(t>0).*(t<2);
inpulse = 1.*(n==1)+1.*(n==3)+1.*(n==5);%1.*(n==5)
figure(1);
subplot(3,1,1);
plot(n,s);

figure(1);
subplot(3,1,2);
plot(n,inpulse);

[conv_y,out_length] = conv_(s,inpulse);
conv_n = linspace(2*min(n),2*max(n),out_length);
figure(1);
subplot(3,1,3);
plot(conv_n,conv_y);

%% LAB 03 Practice 1 Q2
%s = cos(2*pi*0.2*t).*(t<-pi) + cos(2*pi*0.8*t).*(t>=-1*pi);
%s = cos(2*pi*0.2*t)+cos(2*pi*0.8*t)+cos(2*pi*0.4*t);
s = cos(2*pi*0.2*n);
zero_syc = (length(n)-1)/2 +1;
h=zeros(1,length(n));
for i = 2:length(n)
    h(i)=h(i-1)+2*(n(i)==1)+3*(n(i)==4);    
end
figure(2);
subplot(3,1,1);
plot(n,s);

figure(2);
subplot(3,1,2);
plot(n,h);

[conv_y,out_length] = conv_(s,h);
conv_n = linspace(2*min(n),2*max(n),out_length);
figure(2);
subplot(3,1,3);
plot(conv_n,conv_y);

%z =  linspace(0,5,1001);
%H_of_z = (1 + 0.9*z.^-1 )./((1-1.6*z.^-1+0.73*z.^-2));
%H_of_z = ((1 - 1*z.^-1 ).*(1 - 4*z.^-1 ))./(((1 - 2*z.^-1 ).*(1 - 2.1*z.^-1 ).*(1 - 2.2*z.^-1 ).*(1 - 2.3*z.^-1 )));
%H_of_z(isnan(H_of_z))=0;

%{
figure(2);
subplot(2,1,1);
plot(z,abs(H_of_z));
title('|H(z)|')
subplot(2,1,2);
plot(z,phase(H_of_z));
title('∠H(z)')

S = fftshift(fft(s));
h_of_n = ifft(H_of_z);

figure(3);
plot(n,abs(h_of_n));

figure(4);
plot(z,abs(S));
%}



