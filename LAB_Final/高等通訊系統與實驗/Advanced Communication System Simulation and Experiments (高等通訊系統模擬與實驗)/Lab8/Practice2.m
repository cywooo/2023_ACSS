clear all;
close all;

%% BPSK sequence
N = 32;
data = randi([0 1],1,N);
data(data==0) = -1;
BPSK_seq = data;

%% SRRC pulse shaping in digital domain
% Up-sampling
Ld = 4;
upsample_digital = zeros(1,N*Ld);
upsample_digital(1:Ld:end) = BPSK_seq;

% Filter
trun = 5;
M = Ld;
alpha = 0.3;
filter_digital = conv(upsample_digital,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
filter_digital = filter_digital(delay+1:end-delay);
filter_digital = filter_digital/max(filter_digital);

%% DMA filter
% Up-sampling
La = 4;
upsample_analog = zeros(1,length(filter_digital)*La);
upsample_analog(1:La:end) = filter_digital;

% Filter
trun = 5;
M = Ld*La;
alpha = 0.3;
filter_analog = conv(upsample_analog,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
filter_analog = filter_analog(delay+1:end-delay);
filter_analog = filter_analog/max(filter_analog);

%% Noise
SNR = 10;
signal_power = 1; % for BPSK
noise = signal_power * 10^(-SNR/10);
received_signal = filter_analog + noise;

%% Received signal recover
% Filter
trun = 5;
M = La*Ld;
alpha = 0.3;
received_filter_analog = conv(received_signal,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
received_filter_analog = received_filter_analog(delay+1:end-delay);
received_filter_analog = received_filter_analog/max(received_filter_analog);

% Downsampling
Da = La;
downsample_analog = received_filter_analog(1:Da:end); % M phases

% Filter
trun = 5;
M = Ld;
alpha = 0.3;
received_filter_digital = conv(downsample_analog,SRRC_filter(trun,M,alpha));
delay = (length(SRRC_filter(trun,M,alpha))-1)/2;
received_filter_digital = received_filter_digital(delay+1:end-delay);
received_filter_digital = received_filter_digital/max(received_filter_digital);

% Downsampling
Dd = Ld;
downsample_digital = received_filter_digital(1:Dd:end); % M phases

%% Plot
figure(1);
subplot(2,1,1);
stem(BPSK_seq);
title('Original symbol');
subplot(2,1,2);
stem(downsample_digital);
title('Recovered symbol');

figure(2);
stem(downsample_digital);
hold on
stem(BPSK_seq);
legend('Recovered','Original');

