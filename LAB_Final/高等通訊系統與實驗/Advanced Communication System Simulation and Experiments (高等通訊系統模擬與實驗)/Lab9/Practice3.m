clear all
close all

%% Transmit
% BPSK
N = 64;
data = randi([0 1],1,N);
data(data==0) = -1;

% Upsample
Lu = 16;
upsample_data = zeros(1,N*Lu);
upsample_data(1:Lu:end) = data;

% Hybrid shaping
trun = 5;
M = Lu;
alpha = 0.3;
hybrid_data = conv(upsample_data,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
hybrid_data = hybrid_data(delay+1:end-delay);

% DAC
Ldac = 2;
dac_data = zeros(1,length(hybrid_data)*Ldac);
dac_data(1:Ldac:end) = hybrid_data;

% DMA
Ldma = Lu*Ldac;
dma_data = filter(IIR_DMA,dac_data);
group_delay = 1;
dma_data = dma_data(group_delay+1:end);

% Up-conversion
f = 8*10^6;
R = 1*10^6;
fs = Ldma*R;
fdc = f/fs;
wc = 2*pi*fdc;
t = 0:length(dma_data)-1;
upconvert_data = real(dma_data.*exp(i*wc*t));

%% Receive (Intermediate-frequency demodulation)(IF)
% Downconversion
f_IF = 4*10^6;
w_IF = 2*pi*f_IF/fs;
downconvert_data = upconvert_data.*cos((wc-w_IF)*t);

% DMA
received_dma = filter(IIR_DMA,downconvert_data);
group_delay = 1;
received_dma = received_dma(group_delay+1:end);

% ADC
Ladc = 2;
adc_data = received_dma(1:Ladc:end);

% Downconversion (IF)
w_IF_2 = 2*pi*f_IF/(Lu*R);
k = 0:length(adc_data)-1;
IF_received = adc_data.*exp(-i*w_IF_2*k);

% Hybrid shaping
trun = 5;
M = Lu;
alpha = 0.3;
hybrid_received = conv(IF_received,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
hybrid_received = hybrid_received(delay+1:end-delay);

% Downsample
Ld = Lu;
downsample_data = hybrid_received(1:Ld:end);

% Detection
detect = (downsample_data>0) - (downsample_data<0);

figure(1);
stem(data);
hold on
stem(detect);
