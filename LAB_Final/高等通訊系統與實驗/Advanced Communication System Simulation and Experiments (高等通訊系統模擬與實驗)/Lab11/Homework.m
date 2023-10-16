clear all
close all

% Input Signal
N = 128;
a = zeros(1,N);
a(1,1) = 1;

% Up-sampling
Lu = 16;
upsample_a = zeros(1,length(a)*Lu);
upsample_a(1:Lu:end) = a;

% Hybrid pulse shaping
trun = 5;
M = Lu;
alpha = 1;
hybrid_a = conv(upsample_a,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
hybrid_a = hybrid_a(delay+1:end-delay);

% DAC
Ldac = 2;
dac_a = zeros(1,length(hybrid_a)*Ldac);
dac_a(1:Ldac:end) = hybrid_a;

% Analog filter (DMA)
Ldma = Lu*Ldac;
dma_a = filter(IIR_DMA,dac_a);
group_delay = 1;
dma_a = dma_a(group_delay+1:end);

% Up-conversion
fc = 8*10^6;
R = 1*10^6;
fs = R*Ldma;
fdc = fc/fs;
wc = 2*pi*fdc;
t = 1:length(dma_a);
upconvert_a = real(dma_a.*exp(i*wc*t));

% Channel
de = Ldma-1;
h = [0.3 zeros(1,de) 1 zeros(1,de) 0.3];
y_h = conv(h,upconvert_a);
y_h = y_h(1:end-de);
figure(1);
stem(h);
title('Channel');

% Noise
SNR = 10;
signal_power = mean(abs(a).^2);
noise_power = signal_power*10^(-SNR/10);
noise = (noise_power/2)*(randn(1,length(y_h))+1i*randn(1,length(y_h)));
y = y_h + noise;

% Down-conversion
k = 1:length(y);
downconvert_a = y.*exp(-i*wc*k);

% Analog filter (DMA)
dma_received = filter(IIR_DMA,downconvert_a);
dma_received = dma_received(group_delay+1:end);

% ADC
Ladc = Ldac;
adc_a = dma_received(1:Ladc:end);

% Hybrid pulse shaping
hybrid_received = conv(adc_a,SRRC_filter(trun,M,alpha));
hybrid_received = hybrid_received(delay+1:end-delay);

% Down-sampling
Ld = Lu;
downsample_a = hybrid_received(1:Ld:end);
figure(2);
stem(real(downsample_a));
hold on
stem(imag(downsample_a));
title('Downsampled Signal');
legend('Real part','Imaginary part');

% Equalization
n = Ldma;
ZF_equalizer_causal = (10/3)*(-1/3).^(0:n-1);
ZF_equalizer_noncausal = -(-3).^(-n:-1);
recover_tmp = conv(downsample_a,ZF_equalizer_causal);
recover = conv(recover_tmp,ZF_equalizer_noncausal);
response_delay = n;
recover = recover(response_delay+1:end-de);

% Plot
figure(3);
subplot(2,1,1);
stem(a);
title('Transmitted Signal');
subplot(2,1,2);
stem(real(recover));
hold on
stem(imag(recover));
title('Recovered Signal');
legend('Real part','Imaginary part');

