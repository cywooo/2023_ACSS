clear all
close all

nt = 2;
nr = 2;
%% BPSK
N = 64;
symbol = randi([0 1],nt,N);
symbol(symbol==0) = -1;

%% Upsample
Lu = 16;
upsample_sym = zeros(nt,length(symbol)*Lu);
upsample_sym(1,1:Lu:end) = symbol(1,:);
upsample_sym(2,1:Lu:end) = symbol(2,:);

%% SRRC
trun = 5;
M = Lu;
alpha = 0.3;
srrc_sym_1 = conv(upsample_sym(1,:),SRRC_filter(trun,M,alpha));
srrc_sym_2 = conv(upsample_sym(2,:),SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym_1 = srrc_sym_1(delay+1:end-delay);
srrc_sym_2 = srrc_sym_2(delay+1:end-delay);

%% DAC
Ldac = 4;
dac_sym = zeros(nt,length(srrc_sym_1)*Ldac);
dac_sym(1,1:Ldac:end) = srrc_sym_1;
dac_sym(2,1:Ldac:end) = srrc_sym_2;

%% DMA
dma_sym_1 = filter(DMA,dac_sym(1,:));
dma_sym_2 = filter(DMA,dac_sym(2,:));
group_delay = 3;
dma_sym_1 = dma_sym_1(group_delay+1:end);
dma_sym_2 = dma_sym_2(group_delay+1:end);

%% Upconversion
f = 8*1e6;
R = 1*1e6;
fs = R*Lu*Ldac;
fdc = f/fs;
m = 0:length(dma_sym_1)-1;
upconvert_sym_1 = real(dma_sym_1 .* exp(1i*2*pi*fdc*m));
upconvert_sym_2 = real(dma_sym_2 .* exp(1i*2*pi*fdc*m));
upconvert_sym = [upconvert_sym_1; upconvert_sym_2];

%% Channel
H = randn(nr,nt);
H = H ./ sqrt(sum(abs(H).^2));
channel_sym = H*upconvert_sym;

%% Noise
SNR = 10;
signal_power_1 = mean(abs(upconvert_sym_1).^2);
noise_power_1 = signal_power_1*10^(-SNR/10);
noise_1 = sqrt(noise_power_1/2) * (randn(1,length(channel_sym))+1i*randn(1,length(channel_sym)));
received_signal_1 = channel_sym(1,:) + noise_1;

signal_power_2 = mean(abs(upconvert_sym_2).^2);
noise_power_2 = signal_power_2*10^(-SNR/10);
noise_2 = sqrt(noise_power_2/2) * (randn(1,length(channel_sym))+1i*randn(1,length(channel_sym)));
received_signal_2 = channel_sym(2,:) + noise_2;

received_signal = [received_signal_1; received_signal_2];

%% Downconversion
m = 0:length(received_signal)-1;
downconvert_sym = received_signal .* exp(1i*2*pi*fdc*m);

%% DMA
dma_received_1 = filter(DMA,downconvert_sym(1,:));
dma_received_2 = filter(DMA,downconvert_sym(2,:));
dma_received_1 = dma_received_1(group_delay+1:end);
dma_received_2 = dma_received_2(group_delay+1:end);

%% ADC
Ladc = Ldac;
adc_received_1 = dma_received_1(1:Ladc:end);
adc_received_2 = dma_received_2(1:Ladc:end);

%% SRRC
srrc_received_1 = conv(adc_received_1,SRRC_filter(trun,M,alpha));
srrc_received_2 = conv(adc_received_2,SRRC_filter(trun,M,alpha));
srrc_received_1 = srrc_received_1(delay+1:end-delay);
srrc_received_2 = srrc_received_2(delay+1:end-delay);

%% Downsample
Ld = Lu;
downsample_received_1 = srrc_received_1(1:Ld:end);
downsample_received_2 = srrc_received_2(1:Ld:end);
downsample_received = [downsample_received_1; downsample_received_2];

%% ZF Detection
ZF_detector = pinv(H);
ZF_recover = real(ZF_detector*downsample_received);

%% Detection
detection = (ZF_recover>0)-(ZF_recover<0);

%% Plot
figure(1);
subplot(2,1,1);
stem(symbol(1,:));
hold on
stem(ZF_recover(1,:));
title('Recover of first signal');
legend('Original','Recovered');
subplot(2,1,2);
stem(symbol(2,:));
hold on
stem(ZF_recover(2,:));
title('Recover of second signal');
legend('Original','Recovered');

figure(2);
subplot(2,1,1);
stem(symbol(1,:));
hold on
stem(detection(1,:));
title('Detection of first signal');
legend('Original','Recovered');
subplot(2,1,2);
stem(symbol(2,:));
hold on
stem(detection(2,:));
title('Detection of second signal');
legend('Original','Recovered');
