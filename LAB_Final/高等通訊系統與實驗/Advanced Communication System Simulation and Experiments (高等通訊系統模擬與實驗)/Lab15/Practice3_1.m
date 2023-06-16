clear all
close all

%% Input signal
N = 1024;
input_bit = randi([0 1],1,N);
symbol = input_bit;
symbol(symbol==0) = -1;

%% Upsampling
Lu = 16;
upsample_sym = zeros(1,length(symbol)*Lu);
upsample_sym(1:Lu:end) = symbol;

%% Gaussian filter
M = Lu;
BT = 0.5;
trun = 5;
gaussian_sym = conv(upsample_sym,Gaussian_filter(M,BT,trun),'same');

%% Summation to phase (Integration)
Sf = zeros(1,length(gaussian_sym));
for m=1:length(gaussian_sym)
    Sf(m) = sum(gaussian_sym(1:m));
end
fd = 150*1e3;
fb = 1*1e6;
Tb = 1/(fb*Lu);
integrate_sym = exp(1i*2*pi*fd*Tb.*Sf);

%% IF modulation
f_IF = 2*1e6;
fs = fb*Lu;
m = 0:length(integrate_sym)-1;
IF_sym = real(integrate_sym .* exp(1i*2*pi*f_IF/fs*m));

%% DAC
Ldac = 4;
dac_sym = zeros(1,Ldac*length(IF_sym));
dac_sym(1:Ldac:end) = IF_sym;

%% DMA
dma_sym = filter(DMA_new,dac_sym);
group_delay = 11;
dma_sym = dma_sym(group_delay+1:end);

%% Noise
SNR = 7;
signal_power = mean(abs(dma_sym).^2);
noise_power = signal_power * 10^(-SNR/10);
noise = sqrt(noise_power/2) * (randn(1,length(dma_sym)) + 1i*randn(1,length(dma_sym)));
received_signal = dma_sym + noise;

%% DMA (Receiver)
dma_received = filter(DMA_new,received_signal);
dma_received = dma_received(group_delay+1:end);

%% ADC
Ladc = Ldac;
adc_sym = dma_received(1:Ladc:end);

%% IF demodulation
m = 0:length(adc_sym)-1;
IF_received = adc_sym .* exp(-1i*2*pi*f_IF/fs*m);

%% SRRC filter
trun = 5;
M = Lu;
alpha = 0.4;
srrc_sym = conv(IF_received,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
srrc_sym = srrc_sym(delay+1:end-delay);

%% Phase
phase_sym = unwrap(angle(srrc_sym))/(2*pi*fd*Tb);
phase_received = zeros(1,length(phase_sym));
phase_received(1) = phase_sym(1);
for m=2:length(phase_sym)
    phase_received(m) = phase_sym(m) - phase_sym(m-1);
end

%% Gaussian filter
M = Lu;
gaussian_received = conv(phase_received,Gaussian_filter(M,BT,trun),'same');

%% Downsampling
Ld = Lu;
downsample_received = gaussian_received(1:Ld:end);

%% Detection
detection = (downsample_received>0) - (downsample_received<0);
recover_bit = detection;
recover_bit(recover_bit==-1) = 0;

%% BER
BER = sum(recover_bit ~= input_bit)/length(recover_bit);
disp(['The BER is ',num2str(BER),' dB']);


