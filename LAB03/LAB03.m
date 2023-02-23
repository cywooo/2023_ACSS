clear all; clc ;close all;
%% LAB 03 Practice 1
%01
t =  linspace(-10,10,101);
s = cos(2*pi*0.2*t);
%s = t.*(t>0).*(t<2);
inpulse = 1.*(t==5);
figure(1);
subplot(3,1,1);
plot(t,s);

figure(1);
subplot(3,1,2);
plot(t,inpulse);

y = conv_(s,inpulse);
n = (-length(y)+1)/2:(length(y)-1)/2;
figure(1);
subplot(3,1,3);
plot(n,y);

%02
%s = cos(2*pi*0.2*t).*(t<-pi) + cos(2*pi*0.8*t).*(t>=-1*pi);
%s = cos(2*pi*0.2*t)+cos(2*pi*0.8*t)+cos(2*pi*0.4*t);
%figure;
%plot(t,s);


