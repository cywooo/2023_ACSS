clear all; clc ;close all;
%{
t = linspace(-1,1,1e3);
f1 = 0.5;
f2 = 1;
f3 = 2;
x1 = cos(2*pi*f1*t);
x2 = cos(2*pi*f2*t);
x3 = cos(2*pi*f3*t);

figure;
hold on
title('Frequency comparing');
plot(t,x1);
plot(t,x2);
plot(t,x3);
hold off
%}

%{
f = 0.01;
x = cos(2*pi*f*[0:100]);
stem(x);
%}

%% Practice 2
power = 4 ;
a = 2 ;
f = a/2^power;
x = cos(2*pi*f*[1:2^power]);
figure(2);
subplot(2,1,1);
stem(x);
title(['Power of 2 is ', num2str(power),' f = ',num2str(a),'/',num2str(2^power)]);
%X = fftshift(fft(x));
X = fft(x);
Z = [0:2^power-1];
subplot(2,1,2);
stem(Z,abs(X));

%% Practice 3
power = 4 ;
a = 0.7 ;
f = a/2^power;
x = cos(2*pi*f*[1:2^power]);
figure(3);
subplot(2,2,1);
stem(x);
title(['Power of 2 is ', num2str(power),' &  f = ',num2str(a),'/',num2str(2^power)]);
%X = fftshift(fft(x));
X = fft(x);
Z = [0:2^power-1];
subplot(2,2,2);
stem(Z,abs(X));


x2 = [x,zeros(1,2^power)];
subplot(2,2,3);
stem(x2);
title(['Power of 2 is ', num2str(power),' &  f = ',num2str(a),'/',num2str(2^power)]);
%X = fftshift(fft(x));
X2 = fft(x2);
Z2 = [0:length(X2)-1];
subplot(2,2,4);
stem(Z2,abs(X2));

a = 1 ;
f = a/2^power;
x = cos(2*pi*f*[1:2^power]);
figure(4);
subplot(2,2,1);
stem(x);
title(['Power of 2 is ', num2str(power),' &  f = ',num2str(a),'/',num2str(2^power)]);
%X = fftshift(fft(x));
X = fft(x);
Z = [0:2^power-1];
subplot(2,2,2);
stem(Z,abs(X));


x2 = [x,zeros(1,2^power)];
subplot(2,2,3);
stem(x2);
title(['Power of 2 is ', num2str(power),' &  f = ',num2str(a),'/',num2str(2^power)]);
%X = fftshift(fft(x));
X2 = fft(x2);
Z2 = [0:length(X2)-1];
subplot(2,2,4);
stem(Z2,abs(X2));

%% Practice 4
a = 5 ;
M = 2^a ;
x4_ = zeros(1,M);
x4_(M/4+1:M*3/4)=1;
x4 = fftshift(x4_);
X4 = fftshift(ifft(x4));
X4_press = zeros()
Z4 = [-length(X4)/2:length(X4)/2-1];
hold on ;
figure(4)
subplot(2,1,1);
stem(Z4,abs(x4_));
title("x");
subplot(2,1,2);
plot(Z4,abs(X4));
title("X");

%% Practice 5
signal_power = 10;
SNR = 10 ;
x = [zeros(1,6),[0:8],7-[0:7],zeros(1,6)];
x_power = mean(x.^2);
x = x./(sqrt(signal_power/x_power));

n = randn(1,length(x))*sqrt(signal_power/(10^(SNR/10)));
s = x + n;
z = [-14:14];
figure(5);
subplot(3,1,1);
plot(z,n);
title('noise');
subplot(3,1,2);
plot(z,x);
title('signal');
subplot(3,1,3);
plot(z,s);
%title(['signal with noise  SNR=',num2str(SNR)]);

n_power = mean(n.^2);
s_power = mean(s.^2);

ww= 5;
n_ = zeros(1,length(n));
x_ = zeros(1,length(n));
n_(15-ww:15+ww)=n(15-ww:15+ww);
x_(15-ww:15+ww)=x(15-ww:15+ww);
figure;
subplot(2,1,1);
plot(z,x_);
title('window signal');
subplot(2,1,2);
plot(z,n_);
title('window noise');
x_w_power = mean(x_.^2);
n_w_power = mean(n_.^2);

10*log10(s_power/n_power)
10*log10(x_w_power/n_w_power)


%{
S = fftshift(fft(s));
figure(6);
subplot(2,1,1);
plot(z,abs(S));

ww= 5;
S_ = zeros(1,length(S));
S_(15-ww:15+ww)=S(15-ww:15+ww);
subplot(2,1,2);
plot(z,abs(S_));

s_ = ifft(fftshift(S_));
figure(7);
plot(z,abs(s_));

e = s_ - s;
figure(8);
plot(z,abs(e));
%}