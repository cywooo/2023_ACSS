clear all; clc; close all;
n = [1:2000];
w = linspace(-1/2,1/2,length(n));

f = 1/8 ;
down = 3;
up = 3;
s = cos(2*pi*f*n);
S = fftshift(fft(s));

figure(1);
subplot(5,1,1);
plot(n,s);
title("s");
figure(2);
subplot(5,1,1);
plot(w,abs(S));
title("|S|");

down_s = s([1:down:length(s)]);

down_S = fftshift(fft(down_s));

figure(1);
subplot(5,1,2);
plot([1:length(down_s)],down_s);
title("s'");
figure(2);
subplot(5,1,2);
plot(linspace(-1/2,1/2,length(down_S)),abs(down_S));
title("|S'|");

up_s = zeros(up,length(down_s)) ;
up_s(1,:) = down_s;
up_s = reshape(up_s,1,up*length(up_s));

up_S = fftshift(fft(up_s));

figure(1);
subplot(5,1,3);
plot([1:length(up_s)],up_s);
title("s'");
figure(2);
subplot(5,1,3);
plot(linspace(-1/2,1/2,length(up_S)),abs(up_S));
title("|S'|");

FIR_up_s = filter(LPF_FIR,up_s);
FIR_up_S = fftshift(fft(FIR_up_s));

figure(1);
subplot(5,1,4);
plot([1:length(FIR_up_s)],FIR_up_s);
title("FIR s'");
figure(2);
subplot(5,1,4);
plot(linspace(-1/2,1/2,length(FIR_up_S)),abs(FIR_up_S));
title("|FIR S'|");

IIR_up_s = filter(LPF_IIR,up_s);
IIR_up_S = fftshift(fft(IIR_up_s));

figure(1);
subplot(5,1,5);
plot([1:length(IIR_up_s)],IIR_up_s);
title("IIR s'");
figure(2);
subplot(5,1,5);
plot(linspace(-1/2,1/2,length(IIR_up_S)),abs(IIR_up_S));
title("|IIR S'|");

MSE_FIR = 10*log10(mean((FIR_up_s([length(FIR_up_s)-length(s):length(FIR_up_s)-1])/0.3-s).^2));
MSE_IIR = 10*log10(mean((-IIR_up_s([length(FIR_up_s)-length(s)+1:length(FIR_up_s)])/0.3-s).^2));

%% DEMO
f_demo = 1/16;
s_demo = cos(2*pi*f_demo*n);
S_demo = fftshift(fft(s_demo));

s_demo_up = up_sample(32,s_demo);
S_demo_up = fftshift(fft(s_demo_up));

s_demo_up_LPF = filter(LPF_FIR_demo,s_demo_up);
S_demo_up_LPF = fftshift(fft(s_demo_up_LPF));

s_demo_down = down_sample(32,s_demo_up_LPF);
S_demo_down = fftshift(fft(s_demo_down));

figure(3);
subplot(4,1,1);
plot([1:length(s_demo)],s_demo);
figure(4);
subplot(4,1,1);
plot(linspace(-1/2,1/2,length(S_demo)),abs(S_demo));
figure(3);
subplot(4,1,2);
plot([1:length(s_demo_up)],s_demo_up);
figure(4);
subplot(4,1,2);
plot(linspace(-1/2,1/2,length(S_demo_up)),abs(S_demo_up));
figure(3);
subplot(4,1,3);
plot([1:length(s_demo_up_LPF)],s_demo_up_LPF);
figure(4);
subplot(4,1,3);
plot(linspace(-1/2,1/2,length(S_demo_up_LPF)),abs(S_demo_up_LPF));
figure(3);
subplot(4,1,4);
plot([1:length(s_demo_down)],s_demo_down);
figure(4);
subplot(4,1,4);
plot(linspace(-1/2,1/2,length(S_demo_down)),abs(S_demo_down));

