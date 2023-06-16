clear all;
close all;
% Conduct QPSK sequence
N = 128;
data = randi([0 1],1,N);
reshape_data = reshape(data,2,N/2);
real_part = ones(1,N/2);
imaginary_part = ones(1,N/2);
real_part(reshape_data(1,:)==0) = -1;
imaginary_part(reshape_data(2,:)==0) = -1;
QPSK_symbol = real_part + 1i*imaginary_part;

% Upsampling
L = 64;
upsampled_sym = zeros(1,L*length(QPSK_symbol));
upsampled_sym(1:L:end) = QPSK_symbol;

% Practical DAC
practical_sym = conv(upsampled_sym,ones(1,L));

% SRRC pulse shaping
truncate_T = 5;
M = 64;
alpha = 0.3;
pulseShaping_sym = conv(practical_sym,SRRC_filter(truncate_T,M,alpha),'same');

% Up-conversion
m = length(practical_sym);
symbol_rate = 1*10^6;
f = 8*10^6;
fs = symbol_rate*L;
f_dc = f/fs;
wc = 2*pi*f_dc;
upconversion_sym = real(pulseShaping_sym .* exp(1i*wc.*[1:m]));
t = (0:length(upconversion_sym)-1)/fs;
figure(1);
plot(t,upconversion_sym);
title('Up-Converted Signal');
freq = fs*((-m/2)+1:(m/2))/m;
figure(2);
plot(freq,fftshift(abs(fft(upconversion_sym))));
title('Up-Converted Spectrum');
