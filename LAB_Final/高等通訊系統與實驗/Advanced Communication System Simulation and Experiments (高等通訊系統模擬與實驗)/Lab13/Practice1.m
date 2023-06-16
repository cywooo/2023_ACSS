clear all
close all

%% QPSK
N = 64;
data = randi([0 1],2,N);
data(data==0) = -1;
symbol = data(1,:) + 1i*data(2,:);

%% Channel
nt = 1;
nr = 2;
h = 1/sqrt(2)*(randn(nr,nt)+1i*randn(nr,nt));
h = h / sqrt(sum(abs(h).^2));
conv_sym = h*symbol;

%% Noise
SNR = 10;
signal_power = mean(abs(conv_sym).^2,'all');
noise_power = signal_power*10^(-SNR/10);
noise = sqrt(noise_power/2)*(randn(nr,length(conv_sym))+1i*randn(nr,length(conv_sym)));
received_signal = conv_sym + noise;

%% Received Beamforming
recover = h'*received_signal;

%% Plot
figure(1);
stem(real(symbol));
hold on
stem(real(recover));
title('Recovered Real Part');
legend('Original','Recovered');

figure(2);
stem(imag(symbol));
hold on
stem(imag(recover));
title('Recovered Imaginary Part');
legend('Original','Recovered');
