clear all
close all

%% Magnitude of the first path is larger than that of the second
%Input signal
N = 128;
a = zeros(1,N);
a(1,1) = 1;

% Up-sampling
Lu = 4;
upsample_a = zeros(1,length(a)*Lu);
upsample_a(1:Lu:end) = a;

% Hybrid shaping
trun = 5;
alpha = 1;
M = Lu;
hyb_a = conv(upsample_a,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
hyb_a = hyb_a(delay+1:end-delay);

% DAC
Ldac = 8;
dac_a = zeros(1,length(hyb_a)*Ldac);
dac_a(1:Ldac:end) = hyb_a;

% DMA
Ldma = Lu*Ldac;
dma_a = filter(IIR_DMA,dac_a);
group_delay = 8;
dma_a = dma_a(group_delay+1:end);

% Up-conversion
fc = 8*10^6;
R = 1*10^6;
fs = R*Ldma;
fdc = fc/fs;
wc = 2*pi*fdc;
t = 0:length(dma_a)-1;
upconvert_a = real(dma_a.*exp(i*wc*t));

% Channel
tao = Ldma-1;
h = [1 zeros(1,tao) -0.5];
h_a = conv(h,upconvert_a);
h_a = h_a(1:end-Ldma+1);
figure(1);
stem(h);
title('Channel');

% % Noise
% SNR = 20;
% signal_power = mean(abs(a).^2);
% noise_power = signal_power*10^(-SNR/10);
% noise = sqrt(noise_power/2)*(randn(1,length(h_a))+1i*randn(1,length(h_a)));
% y = h_a + noise;

% Down-conversion
y = h_a;
t = 0:length(y)-1;
downconvert_a = y.*exp(-i*wc*t);

% DMA
dma_received = filter(IIR_DMA,downconvert_a);
dma_received = dma_received(group_delay+1:end);

% ADC
Ladc = Ldac;
adc_a = dma_received(1:Ladc:end);

% Hybrid shaping
hyb_received = conv(adc_a,SRRC_filter(trun,M,alpha));
hyb_received = hyb_received(delay+1:end-delay);

% Down-sampling
Ld = Lu;
downsample_a = hyb_received(1:Ld:end);

% Plot
figure(2);
subplot(2,1,1);
stem(a);
title('Transmit Signal');
subplot(2,1,2);
stem(real(downsample_a));
hold on
stem(imag(downsample_a));
title('Received Signal');
