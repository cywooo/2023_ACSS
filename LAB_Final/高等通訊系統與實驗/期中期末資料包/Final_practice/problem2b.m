clear all; close all; clc;
data = load('data_given');
x = data.a_2a;
intf = data.intf_2b;

% Tx
upsample1 = 16;
x_u1 = zeros(1,length(x)*upsample1);
x_u1(1:upsample1:end) = x;

alpha = 0.3;
span = 5;

x_SRRC_Tx = conv(x_u1,rcosine(1,16,'fir/sqrt',0.3,5));

upsample2 = 8;
x_u2 = zeros(1,length(x_SRRC_Tx)*upsample2);
x_u2(1:upsample2:end) = x_SRRC_Tx;

figure(1)
fx = linspace(-1/2,1/2,length(fftshift(fft(x_u2))));
plot(fx,abs(fftshift(fft(abs(x_u2)))))

x_DMA_Tx = filter(IIR_2a_1,x_u2);
group_delay = 14;
x_DMA_Tx = x_DMA_Tx(group_delay:end);

% x_DMA_Tx = filter(IIR_2a,x_u2);


figure(2)
fx = linspace(-1/2,1/2,length(abs(fftshift(fft(x_DMA_Tx)))));
plot(fx,fftshift(abs(fft((x_DMA_Tx)))))

R = 1*10^6;
fc = 32*10^6;
M = upsample1*upsample2;
fdc = fc/(R*M);

m = 1:length(x_DMA_Tx); %%%%%%%
x_IF_Tx = x_DMA_Tx.*exp(1j*2*pi*fdc*m);
x_Tx = real(x_IF_Tx);

figure(3)
fx = linspace(-1/2,1/2,length(fftshift(fft(x_Tx))));
plot(fx,abs(fftshift(fft(x_Tx))))

%% Rx
% adding interference 
y_n = x_Tx + intf(1:end-13);

figure(4)
fx = linspace(-1/2,1/2,length(fftshift(abs(fft(y_n)))));
plot(fx,fftshift(abs(fft(y_n))));

y_irf = filter(imag_filter,y_n);
g_d = 55;
y_irf = y_irf(g_d:end);

figure(5)
fx = linspace(-1/2,1/2,length(fftshift(fft(y_irf))));
plot(fx,abs(fftshift(fft(y_irf))))

f_IF = 2*10^6;
m = 1:length(y_irf);
y_IF = y_irf.*cos( 2*pi*(fdc-f_IF/(R*M)) *m);

y_dma = filter(IIR_2b_rx,y_IF);
group_delay = 8;
y_dma = y_dma(group_delay:end);


% y_dma = filter(IIR_2a,y_IF);
% group_delay = 6;
% y_dma = y_dma(group_delay+1:end);


figure(6)
fx = linspace(-1/2,1/2,length(fftshift(fft(y_dma))));
plot(fx,abs(fftshift(fft(y_dma))))

downsample1 = upsample2;
y_down1 = y_dma(1:downsample1:end);
k = 1:length(y_down1); %%%%%%%
y_IF = y_down1.*exp(-1j*2*pi*( f_IF/(R*upsample1) )*k); %%%--->

figure(7)
fx = linspace(-1/2,1/2,length(fftshift(fft(y_IF))));
plot(fx,abs(fftshift(fft(y_IF))))

% y_SRRC = conv(y_IF,rcosdesign(alpha,span,upsample1),'same');
y_SRRC = conv(y_IF,rcosine(1,16,'fir/sqrt',0.3,5),'same');
% srrc_delay = 80*2;
% y_SRRC = y_SRRC(srrc_delay+1:end);

downsample2 = upsample1;

y = y_SRRC(1:downsample2:end);
y = y(1:end-10)
a_2ha = y;
figure(8)
stem(x);
hold on
stem(y)


% SNR
e = y-x;
s_power = mean(abs(x).^2);
e_power = mean(abs(e).^2);
SNR_2a  = 10*log10(s_power/e_power)
