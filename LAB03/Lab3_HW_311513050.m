clear all; clc ;close all;

f1 = 1/4;
f2 = 1/8;

n = 1:201; % Generate n series
w = linspace(-1/2,1/2,length(n)); % Generate w series
s = 1*cos(2*pi*f1*n) + cos(2*pi*f2*n); % Generate combined signal

%% LAB3 HW1-1 IIR
figure(1); % compare different frequency and output waveform 
subplot(3,1,3);
plot(n,s);
title('Output Signal');
subplot(3,1,1);
plot(n,cos(2*pi*1/4*n));
title('Frequency = ',f1);
subplot(3,1,2);
plot(n,cos(2*pi*1/8*n));
title('Frequency = ',f2);

S = fftshift(fft(s)); 
figure(2);
subplot(3,1,1);
plot(w,abs(S));
title('|S(w)|');

% Generate IIR Filter
h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=0*h(max(j-1,1)).*(j>1)-0.8100*h(max(j-2,1)).*(j>2)+...
        1.*(j == 0+1)-2.4142.*(j == 1+1)+2.4142.*(j == 2+1)-1.*(j == 3+1);       
end
%IIR Filter parameter
% zero on f = 1/8 & f = 0
% - poly([(cos(1/8*pi*2)+sin(1/8*pi*2)*i),(cos(1/8*pi*2)-sin(1/8*pi*2)*i),1]) --> ans=  1.0000   -2.4142    2.4142   -1.0000
% pole near f = 1/4
% - poly([(0.9*cos(1/4*pi*2)+0.9*sin(1/4*pi*2)*i),(0.9*cos(1/4*pi*2)-0.9*sin(1/4*pi*2)*i)]) --> ans=  1.0000   -0.0000    0.8100
figure(6);
polarscatter([1/8*pi,-1/8*pi,0],[1,1,1],120,"O");hold on;
polarscatter([1/4*pi,-1/4*pi],[0.9,0.9],120,"X");hold off;
axis([0 360 0 1.2]);
title('IIR Filter');
legend('Zeros','Poles');

H = fftshift(fft(h));
figure(2);
subplot(3,1,2);
plot(linspace(-1/2,1/2,length(H)),abs(H));
title('|H(w)| IIR');
figure(3);
subplot(4,1,1);
plot(n,s);
title('s(n)');
figure(3);
subplot(4,1,2);
plot(n,h);
title('h(n) IIR');
figure(3);
subplot(4,1,4);
plot(n,cos(2*pi*1/4*n));
axis([0 450 -2 2]);
title('Frequency = ',f1);

%signal pass through filter
[conv_y,out_length] = conv_(s,h);
conv_n = linspace(2*min(n),2*max(n),out_length);
figure(3);
subplot(4,1,3);
plot(conv_n,conv_y);
title("s'(n)");

S_filter = fftshift(fft(conv_y([1:length(n)])));
figure(2);
subplot(3,1,3);
plot(w,abs(S_filter));
title("|S'(w)|");

%% IIR SNR 
%take the power in -(f1+f2)/2 < f < (f1+f2)/2 as noise 
decision = (f1+f2)/2;
SNR_IIR_out = 10*log10((sum((abs(S_filter).*(w(find(S_filter))<-decision)).^2)+sum((abs(S_filter).*(w(find(S_filter))>decision).^2)))/...
                        sum((abs(S_filter).*(w(find(S_filter))>=-decision).*(w(find(S_filter))<=decision)).^2));

%% LAB3 HW1-2 FIR

S = fftshift(fft(s));
figure(4);
subplot(3,1,1);
plot(w,abs(S));
title('|S(w)|');

% Generate FIR Filter
h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=1.*(j == 0+1)-2.4142.*(j == 1+1)+2.4142.*(j == 2+1)-1.*(j == 3+1);       
end
%FIR Filter parameter
%zero on f = 1/8
%poly([(cos(1/8*pi*2)+sin(1/8*pi*2)*i),(cos(1/8*pi*2)-sin(1/8*pi*2)*i),1]) --> ans=  1.0000   -2.4142    2.4142   -1.0000
figure(7);
polarscatter([1/8*pi,-1/8*pi,0],[1,1,1],120,"O");
axis([0 360 0 1.2]);
title('FIR Filter');
legend('Zeros');

H = fftshift(fft(h));
figure(4);
subplot(3,1,2);
plot(linspace(-1/2,1/2,length(H)),abs(H));
title('|H(w)| FIR');
figure(5);
subplot(4,1,1);
plot(n,s);
title('s(n)');
figure(5);
subplot(4,1,2);
plot(n,h);
title('h(n) FIR');
figure(5);
subplot(4,1,4);
plot(n,cos(2*pi*1/4*n));
axis([0 450 -2 2]);
title('Frequency = ',f1);

%signal pass through filter
[conv_y,out_length] = conv_(s,h);
conv_n = linspace(2*min(n),2*max(n),out_length);
figure(5);
subplot(4,1,3);
plot(conv_n,conv_y);
title("s'(n)");

S_filter = fftshift(fft(conv_y([1:length(n)])));
figure(4);
subplot(3,1,3);
plot(w,abs(S_filter));
title("|S'(w)|");

%% FIR SNR 
%take the power in -(f1+f2)/2 < f < (f1+f2)/2 as noise 
decision = (f1+f2)/2;
SNR_FIR_out = 10*log10((sum((abs(S_filter).*(w(find(S_filter))<-decision)).^2)+sum((abs(S_filter).*(w(find(S_filter))>decision).^2)))/...
                        sum((abs(S_filter).*(w(find(S_filter))>=-decision).*(w(find(S_filter))<=decision)).^2));
fprintf('Final SNR output :\n');
fprintf('FIR_SNR = %3.2f , IIR_SNR = %3.2f \n\n',SNR_FIR_out,SNR_IIR_out);