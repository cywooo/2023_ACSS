clear all
close all

%% BPSK
N = 64;
symbol = randi([0 1],1,N);
symbol(symbol==0) = -1;

%% Up-sample
Lu = 16;
upsample_sym = zeros(1,length(symbol)*Lu);
upsample_sym(1:Lu:end) = symbol;

%% SRRC
trun = 5;
M = Lu;
alpha = 0.3;
srrc_sym = conv(upsample_sym,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym = srrc_sym(delay+1:end-delay);

%% DAC
Ldac = 4;
dac_sym = zeros(1,length(srrc_sym)*Ldac);
dac_sym(1:Ldac:end) = srrc_sym;

%% DMA
Ldma = Lu*Ldac;
dma_sym = filter(DMA,dac_sym);
group_delay = 3;
dma_sym = dma_sym(group_delay+1:end);

%% Up-conversion
f = 8*1e6;
R = 1*1e6;
fs = R*Ldma;
fdc = f/fs;
m = 0:length(dma_sym)-1;
upconvert_sym = real(dma_sym .* exp(1i*2*pi*fdc*m));

%% Channel
nt = 1;
nr = 2;
H = randn(nr,nt);%1/sqrt(2) * (randn(nr,nt) + 1i*randn(nr,nt));
H = H/sqrt(sum(abs(H).^2));
channel_sym = H * upconvert_sym;

%% Noise
SNR = 10;
signal_power_1 = mean(abs(upconvert_sym).^2);
noise_power_1 = signal_power_1*10^(-SNR/10);
noise1 = sqrt(noise_power_1/2) * (randn(1,length(channel_sym(1,:)))+1i*randn(1,length(channel_sym(1,:))));
received_signal1 = channel_sym(1,:) + noise1;

signal_power_2 = mean(abs(upconvert_sym).^2);
noise_power_2 = signal_power_2*10^(-SNR/10);
noise2 = sqrt(noise_power_2/2) * (randn(1,length(channel_sym(2,:)))+1i*randn(1,length(channel_sym(2,:))));
received_signal2 = channel_sym(2,:) + noise2;

received_signal = [received_signal1;received_signal2];

%% Down-conversion
m = 0:length(received_signal)-1;
downconvert_sym = received_signal .* exp(-1i*2*pi*fdc*m);

%% DMA
dma_received_1 = filter(DMA,downconvert_sym(1,:));
dma_received_1 = dma_received_1(group_delay+1:end);

dma_received_2 = filter(DMA,downconvert_sym(2,:));
dma_received_2 = dma_received_2(group_delay+1:end);

%% ADC
Ladc = Ldac;
adc_received_1 = dma_received_1(1:Ladc:end);
adc_received_2 = dma_received_2(1:Ladc:end);

%% SRRC
srrc_received_1 = conv(adc_received_1,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_received_1 = srrc_received_1(delay+1:end-delay);

srrc_received_2 = conv(adc_received_2,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_received_2 = srrc_received_2(delay+1:end-delay);

%% Downsample
Ld = Lu;
downsample_received_1 = srrc_received_1(1:Ld:end);
downsample_received_2 = srrc_received_2(1:Ld:end);
downsample_received = [downsample_received_1;downsample_received_2];

%% Received Beamforming
recovered_signal = H'*downsample_received;

%% Plot
figure(1);
stem(symbol);
hold on
stem(real(recovered_signal));

