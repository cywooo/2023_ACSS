clear all
close all

%% QPSK signal
% Generate signals
N = 128;
sI = randi([0 1],1,N);
sI(sI==0) = -1;
sQ = randi([0 1],1,N);
sQ(sQ==0) = -1;
s = sI + 1i*sQ;

% Upsample
Lu = 16;
upsample_s = zeros(1,length(s)*Lu);
upsample_s(1:Lu:end) = s;

% Hybrid pulse shaping
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
gd = 11;
dma_s = dma_s(gd+1:end);

% IQ 
aI = real(dma_s);
aQ = imag(dma_s);
fc = 16*10^6;
R = 1*10^6;
fs = Ldma*R;
fdc = fc/fs;
wc = 2*pi*fdc;
%g = 1.5;
%phi = 60/180*pi;
g = 1;
phi = 0;
m = 0:length(dma_s)-1;
xI = aI*sqrt(2).*cos(wc*m);
xQ = aQ*(-1)*g*sqrt(2).*sin(wc*m+phi);
x = xI + xQ;

% Receiver
y = x;
fi = 4*10^6;
w_IF = 2*pi*fi/fs;
y_IF = y.*cos((wc-w_IF)*m);

% DMA
y_dma = filter(IIR_DMA,y_IF);
y_dma = y_dma(gd+1:end);

% ADC
Ladc = Ldac;
adc_s = y_dma(1:Ladc:end);

% IF demodulation
w_IF_2 = 2*pi*fi/(Lu*R);
k = 0:length(adc_s)-1;
IF_I = adc_s*sqrt(2).*cos(w_IF_2*k);
IF_Q = adc_s*(-1)*sqrt(2).*sin(w_IF_2*k);

% Hybrid pulse shaping
hps_I = conv(IF_I,SRRC_filter(trun,L,alpha));
hps_Q = conv(IF_Q,SRRC_filter(trun,L,alpha));
delay = (length(SRRC_filter(trun,L,alpha))-1)/2;
hps_I = hps_I(delay+1:end-delay);
hps_Q = hps_Q(delay+1:end-delay);
hps_receive = hps_I + 1i*hps_Q;

% Downsample
Ld = Lu;
downsample_s = hps_receive(1:Ld:end);

% Detection
real_detect = (real(downsample_s)>0)-(real(downsample_s)<0);
imag_detect = (imag(downsample_s)>0)-(imag(downsample_s)<0);

figure(1);
stem(sI);
hold on
stem(real_detect);
title('Detection of real part');
legend('Origin','Recovered');
figure(2);
stem(sQ);
hold on
stem(imag_detect);
title('Detection of imaginary part');
legend('Origin','Recovered');

% Verify
alp = (1+g*exp(i*phi))/2;
beta = (1-g*exp(i*phi))/2;
veri_s = (1/2)*(alp*s + beta*conj(s));
figure(3);
plot(downsample_s*4,'o');
hold on
plot(veri_s,'o');
title('Constellation');
legend('Origin','Verify');