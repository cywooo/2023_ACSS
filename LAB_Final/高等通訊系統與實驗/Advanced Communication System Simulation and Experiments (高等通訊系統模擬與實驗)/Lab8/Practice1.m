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
filter_digital = filter_digital(delay:end-delay);

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
filter_analog = filter_analog(delay:end-delay);

figure(1);
subplot(5,1,1);
fs = 1;
f = 2*fs*(-N/2:N/2-1)/N;
plot(f,fftshift(abs(fft(BPSK_seq))));
title('Origin fft');
subplot(5,1,2);
m = length(upsample_digital);
f = 2*fs*(-m/2:m/2-1)/m;
plot(f,fftshift(abs(fft(upsample_digital))));
title('Upsample_1 fft');
subplot(5,1,3);
m = length(filter_digital);
f = 2*fs*(-m/2:m/2-1)/m;
plot(f,fftshift(abs(fft(filter_digital))));
title('filter_1 fft');
subplot(5,1,4);
m = length(upsample_analog);
f = 2*fs*(-m/2:m/2-1)/m;
plot(f,fftshift(abs(fft(upsample_analog))));
title('Upsample_2 fft');
subplot(5,1,5);
m = length(filter_analog);
f = 2*fs*(-m/2:m/2-1)/m;
plot(f,fftshift(abs(fft(filter_analog))));
title('filter_2 fft');
