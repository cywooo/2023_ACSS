clear all
close all

%% Transmitter
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
Ldac = 4;
dac_data = zeros(1,length(hybrid_data)*Ldac);
dac_data(1:Ldac:end) = hybrid_data;

% DMA
Ldma = Lu*Ldac;
dma_data = filter(IIR_DMA_hw,dac_data);
group_delay = 10;
dma_data = dma_data(group_delay+1:end);
figure(1);
ld = length(dma_data);
f_axis = (-ld/2:ld/2-1)/ld;
plot(f_axis,fftshift(abs(fft(dma_data))));
title('Spectrum after DMA Signal');

% Up-conversion
fc = 16*10^6;
R = 1*10^6;
fs = Ldma*R;
fdc = fc/fs;
wc = 2*pi*fdc;
t = 0:length(dma_data)-1;
upconvert_data = real(dma_data.*exp(i*wc*t));
figure(2);
lt = length(t);
f_axis = fs*(-lt/2:lt/2-1)/lt;
plot(f_axis,fftshift(abs(fft(upconvert_data))));
title('Spectrum of Up-converted Signal');

%% Receiver (Intermediate-frequency demodulation)(IF)(Without noise)
% Add Interference on image
fi = 4*10^6;
n = 1:length(upconvert_data);
interference = cos(2*pi*((2*fi-fc)/fs).*n);
received_signal = upconvert_data + interference;
figure(3);
ln = length(n);
f_axis = fs*(-ln/2:ln/2-1)/ln;
plot(f_axis,fftshift(abs(fft(received_signal))));
title('Spectrum of Up-converted Signal adding Interference');

% Image Rejection Filter
image_rejection_data = filter(image_rejection_filter,received_signal);
group_delay = 16;
image_rejection_data = image_rejection_data(group_delay+1:end);
figure(4);
it = length(image_rejection_data);
f_axis = fs*(-it/2:it/2-1)/it;
plot(f_axis,fftshift(abs(fft(image_rejection_data))));
title('Spectrum of Signal passing Image Rejection filter without noise');

% Downconversion
fi = 4*10^6;
w_IF = 2*pi*fi/fs;
m = 0:length(image_rejection_data)-1;
downconvert_data = image_rejection_data.*cos((wc-w_IF)*m);
figure(5);
id = length(downconvert_data);
f_axis = fs*(-id/2:id/2-1)/id;
plot(f_axis,fftshift(abs(fft(downconvert_data))));
title('Spectrum of Downconvert Signal without noise');

% DMA
received_dma = filter(IIR_DMA_hw,downconvert_data);
group_delay = 10;
received_dma = received_dma(group_delay+1:end);

% ADC
Ladc = 4;
adc_data = received_dma(1:Ladc:end);

% Downconversion (IF)
w_IF_2 = 2*pi*fi/(Lu*R);
k = 0:length(adc_data)-1;
IF_received = adc_data.*exp(-i*w_IF_2*k);
figure(6);
ik = length(IF_received);
f_axis = (-ik/2:ik/2-1)/ik;
plot(f_axis,fftshift(abs(fft(IF_received))));
title('Spectrum of IF Downconvert Signal without noise');

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

figure(7);
stem(data);
hold on
stem(detect);
title('Comparison of Recovered signal and Origin signal without adding noise');
legend('Origin','Recovered');

%% Receiver (Intermediate-frequency demodulation)(IF)(With noise)
% Add Interference on image
fi = 4*10^6;
n = 1:length(upconvert_data);
interference = cos(2*pi*((2*fi-fc)/fs).*n);
received_signal_tmp = upconvert_data + interference;

% Add noise
SNR = 10;
signal_power = mean(abs(data).^2);
noise_power = signal_power*10^(-SNR/10);
noise = sqrt(noise_power/2)*(randn(1,length(received_signal_tmp))+1i*randn(1,length(received_signal_tmp)));
received_signal = received_signal_tmp + noise;
figure(8);
rt = length(received_signal);
f_axis = fs*(-rt/2:rt/2-1)/rt;
plot(f_axis,fftshift(abs(fft(received_signal))));
title('Spectrum of Received Signal with noise');

% Image Rejection Filter
image_rejection_data = filter(image_rejection_filter,received_signal);
group_delay = 16;
image_rejection_data = image_rejection_data(group_delay+1:end);
figure(9);
it = length(image_rejection_data);
f_axis = fs*(-it/2:it/2-1)/it;
plot(f_axis,fftshift(abs(fft(image_rejection_data))));
title('Spectrum of Signal passing Image Rejection filter with noise');

% Downconversion
fi = 4*10^6;
w_IF = 2*pi*fi/fs;
m = 0:length(image_rejection_data)-1;
downconvert_data = image_rejection_data.*cos((wc-w_IF)*m);

% DMA
received_dma = filter(IIR_DMA_hw,downconvert_data);
group_delay = 10;
received_dma = received_dma(group_delay+1:end);

% ADC
Ladc = 4;
adc_data = received_dma(1:Ladc:end);

% Downconversion (IF)
w_IF_2 = 2*pi*fi/(Lu*R);
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

figure(10);
stem(data);
hold on
stem(detect);
title('Comparison of Recovered signal and Origin signal adding noise');
legend('Origin','Recovered');
