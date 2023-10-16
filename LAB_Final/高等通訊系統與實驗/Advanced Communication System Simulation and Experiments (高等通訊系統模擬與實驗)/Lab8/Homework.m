clear all
close all

%% Conduct the transmit and receive operation
% BPSK
N = 32;
data = randi([0 1],1,N);
data(data==0) = -1;
BPSK_sym = data;

% Up-sampling (transmit)
Ld = 8;
upsample_digital = zeros(1,N*Ld);
upsample_digital(1:Ld:end) = BPSK_sym;

% Digital pulse shaping (transmit)
trun = 5;
M = Ld;
alpha = 0.3;
digit_ps = conv(upsample_digital,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
digit_ps = digit_ps(delay+1:end-delay);
digit_ps = digit_ps/max(digit_ps);

% Practical DAC (transmit)
La = 4;
upsample_analog_tmp = zeros(1,length(digit_ps)*La);
upsample_analog_tmp(1:La:end) = digit_ps;
upsample_analog = conv(upsample_analog_tmp,ones(1,La));
upsample_analog = upsample_analog(1:length(upsample_analog_tmp));

% DMA filter (transmit)
iir = IIR_DMA_hw;
sos = iir.SOS;
analog_dma = sosfilt(sos,upsample_analog);
%analog_dma = filter(IIR_DMA_hw,upsample_analog);
delay = 20;
analog_dma = analog_dma(delay+1:end);
analog_dma = analog_dma/max(analog_dma);

% Add noise
% SNR = 15;
% signal_power = 1; % for BPSK
% noise_power = signal_power*10^(-SNR/10);
% noise = sqrt(noise_power/2)*(randn(1,length(analog_dma)) + 1i*randn(1,length(analog_dma)));
% received_signal = analog_dma + noise;
received_signal = analog_dma;

% DMA filter (received)
iir = IIR_DMA_hw;
sos = iir.SOS;
received_dma = sosfilt(sos,received_signal);
%received_dma = filter(IIR_DMA_hw,received_signal);
delay = 20;
received_dma = received_dma(delay+1:end);
received_dma = received_dma/max(received_dma);

% ADC (received)
Da = La;
downsample_adc = received_dma(1:Da:end);

% Digital pulse shaping (received)
M = Ld;
rec_ps = conv(downsample_adc,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
rec_ps = rec_ps(delay+1:end-delay);
rec_ps = rec_ps/max(rec_ps);

% Downsampling (received)
Dd = Ld;
recovered_sym = rec_ps(1:Dd:end);
figure(1);
stem(real(recovered_sym));
hold on
stem(BPSK_sym);
legend('Recovered','Origin');

%% Conduct the up-conversion operation
fc = 8*10^6;
L = La*Ld;
R = 1*10^6;
fdc = fc/L/R;
wc = 2*pi*fdc;
fs = L*R;
m = 1:length(analog_dma);
upconvert_sym = real(analog_dma .* exp(1i*wc.*m));
figure(2);
subplot(2,1,1);
n = length(upconvert_sym);
t = (0:n-1)/fs;
plot(t,upconvert_sym);
title('Up-conversion Signal');
subplot(2,1,2);
f = fs*(-n/2:n/2-1)/n;
plot(f,fftshift(abs(fft(upconvert_sym))));
title('Up-conversion Spectrum');
