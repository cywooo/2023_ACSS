% n = 512;
% sd=1;
% t=[1:n];
% fs = 10;
% fq_1 = 1/64;
% fq_2 = 1/4;
% lpf=ones(1,fs)/fs;% b=[0.33 0.33 0.33]
% x = cos(2*pi*fq_1*t);
% ca = sqrt(2)*cos(2*pi*fq_2*t);
% w=sd*randn(1,n);
% mx=x.*ca;   %mod
% mxn=mx+w;
% cmx=mxn.*ca;  %demo
% figure(1)
% zplane(lpf,1);
% figure(2)
% dmx=filter(lpf,1,cmx);
% plot([x' dmx'])






t=[1:512];
fs = 1/64;
fc = 1/4;
%lpf=ones(1,fs)/fs;% b=[0.33 0.33 0.33]
windowsize = 10;
b = (1/windowsize)*ones(1,windowsize);
% b=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];
x = cos(2*pi*fs*t);
ca = sqrt(2)*cos(2*pi*fc*t);

mx=x.*ca;   %mod
cmx=mx.*ca;  %demo

figure(1)
zplane(b,1);
figure(2)
freqz(b,1);
figure(3)
dmx=filter(b,1,cmx);
plot([x' dmx'])




