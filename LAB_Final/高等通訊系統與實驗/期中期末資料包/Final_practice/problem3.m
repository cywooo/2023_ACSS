clear all; close all; clc;
data = load('data_given');
x = data.x_3;
a = data.a_3b;

%% (a)
R = 1*10^6;
L1 = 16; L2 = 8;
fc = 32*10^6;
fdc = fc/(R*L1*L2);
wdc = 2*pi*fdc;

m = 1:length(x);
y_IF = x.*exp(-1j*wdc*m);

figure(1);
f_y_IF = fftshift(fft(y_IF));
fx = linspace(-1/2,1/2,length(f_y_IF));
plot(fx,abs(f_y_IF))

% y_DMA = filter(IIR_3,y_IF);
% g_d = 3;
% y_DMA = y_DMA(g_d:end);
y_DMA = filter(IIR_3a,y_IF);
g_d = 9;
y_DMA = y_DMA(g_d:end);

figure(2);
f_y_DMA = fftshift(fft(y_DMA));
fx = linspace(-1/2,1/2,length(f_y_DMA));
plot(fx,abs(f_y_DMA))

D1 = L2;
y_down1 = y_DMA(1:D1:end);

y_SRRC = conv(y_down1,rcosine(1,16,'fir/sqrt',0.3,5));
%%%%%%%%%%
delay  = 80*2;
y_SRRC = y_SRRC(delay+1:end);
%%%%%%%%%%%

D2 = L1;
y = y_SRRC(1:D2:end);
%%%%%%%%%%
y = y(1:end-10);
%%%%%%%%%%%%
a_3ha = y;

figure(3);
plot(y,'o')

figure(4);
f_y = fftshift(fft(y));
fx = linspace(-1/2,1/2,length(f_y));
plot(fx,abs(f_y))

%% (b)
y_train = y(1:end);
a_train = a(1:end);
y_angle = unwrap(angle(y_train));
a_angle = unwrap(angle(a_train));

e = y_angle-a_angle;
CFO_3b = mean(e)

figure(5)
plot(e)

%% (c)
fc = 32*10^6;
fdc = fc/(R*L1*L2);
wdc = 2*pi*fdc;

m = 1:length(x);
y_IF_com = x.*exp(-1j*(wdc*m-CFO_3b));

y_DMA_com = filter(IIR_3a,y_IF_com);
g_d = 9;
y_DMA_com = y_DMA_com(g_d:end);

D1 = L2;
y_down1_com = y_DMA_com(1:D1:end);

y_SRRC_com = conv(y_down1_com,rcosine(1,16,'fir/sqrt',0.3,5));
%%%%%%%%%%
delay  = 80*2;
y_SRRC_com = y_SRRC_com(delay+1:end);
%%%%%%%%%%%

D2 = L1;
y_com = y_SRRC_com(1:D2:end);
%%%%%%%%%%
y_com = y_com(1:end-10);
%%%%%%%%%%%%
a_3hc = y_com;

figure(6);
plot(y_com,'o')