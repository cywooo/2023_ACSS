clear all
close all
%% QPSK signal
% Generate signals
N = 256;
sI = randi([0 1],1,N);
sI(sI==0) = -1;
sQ = randi([0 1],1,N);
sQ(sQ==0) = -1;
s = sI + i*sQ;

% Upsample
Lu = 16;
upsample_s = zeros(1,length(s)*Lu);
upsample_s(1:Lu:end) = s;

% hybrid pulse shaping
trun = 5;
L = 16;
alpha = 0.3;
hps_s = conv(upsample_s,SRRC_filter(trun,L,alpha));
delay = (length(SRRC_filter(trun,L,alpha))-1)/2;
hps_s = hps_s(delay+1:end-delay);

% DAC
Ldac = 4;
dac_s = zeros(1,length(hps_s)*Ldac);
dac_s(1:Ldac:end) = hps_s;

% DMA
Ldma = 64;
dma_s = filter(IIR_DMA,dac_s);
gd = 10;
dma_s = dma_s(gd+1:end);

% IQ 
aI = real(dma_s);
aQ = imag(dma_s);
fc = 16*10^6;
R = 1*10^6;
fs = Ldma*R;
fdc = fc/fs;
g = 1;
phi = 0;
m = 0:length(dma_s)-1;
xI = aI*sqrt(2).*cos(2*pi*fdc*m);
xQ = aQ*(-1)*g*sqrt(2).*sin(2*pi*fdc*m+phi);
x = xI + xQ;

% CFO
offset = 1000/fs;
mI = x*sqrt(2).*cos(2*pi*(fdc+offset)*m);
mQ = x*(-1)*sqrt(2).*sin(2*pi*(fdc+offset)*m);

% DMA
yI = filter(IIR_DMA,mI);
yI = yI(gd+1:end);
yQ = filter(IIR_DMA,mQ);
yQ = yQ(gd+1:end);
y = yI + i*yQ;

% ADC
Ladc = Ldac;
adc_s = y(1:Ladc:end);

% hybrid pulse shaping
hps_receive = conv(adc_s,SRRC_filter(trun,L,alpha));
delay = (length(SRRC_filter(trun,L,alpha))-1)/2;
hps_receive = hps_receive(delay+1:end-delay);

% Downsample
Ld = Lu;
downsample_s = hps_receive(1:Ld:end);

figure(1);
stem(sI);
hold on
stem(real(downsample_s));
figure(2);
stem(sQ);
hold on
stem(imag(downsample_s));

figure(3);
plot(downsample_s,'o');
hold on
plot(s,'o');