clear all;
close all;
%% Generate a BPSK sequence
N = 64;
data = randi([0 1],1,N);
symbol = ones(1,N);
symbol(data==0) = -1;
figure(1);
stem(symbol);
title('Input BPSK Symbol');

%% Conduct the RC pulse shaping operation with an ideal DAC (upsampled by a factor 32)
% Upsampling
L = 32;
upsample_sym = zeros(1,N*L);
upsample_sym(1:L:end) = symbol;

% Pulse shaping(RC)
alpha = 0.3;
M = 8;
truncate_T = 5;
pulseShaping_s = conv(upsample_sym,RC_filter(truncate_T,M,alpha)); % Call RC_filter function
figure(2);
stem(pulseShaping_s);
title('Pulse Shaping Symbol for ideal DAC');

% Downsampling
D = 32;
downsample_s = pulseShaping_s(1:D:end);
figure(3);
stem(downsample_s);
title('Downsampled Symbol for ideal DAC');

%% Conduct the SRRC pulse shaping operation with an practical DAC (upsampled by a factor 32)
% Upsampling
L = 32;
upsample_sym = zeros(1,N*L);
upsample_sym(1:L:end) = symbol;

% Rect DMA
rect_DMA = conv(upsample_sym,ones(1,L));

% Pulse shaping in TX (SRRC)
alpha = 0.3;
M = 8;
truncate_T = 5;
pulseShaping_s = conv(rect_DMA,SRRC_filter(truncate_T,M,alpha),'same'); % Call SRRC_filter function
figure(4);
stem(pulseShaping_s);
title('Pulse Shaping Symbol for practical DAC');

% Pulse shaping in RX (SRRC)
pulseShaping_s_rx = conv(pulseShaping_s,SRRC_filter(truncate_T,M,alpha),'same'); % Call SRRC_filter function
figure(5);
stem(pulseShaping_s_rx);
title('Received Pulse Shaping Symbol for practical DAC');

% Downsampling
D = 32;
downsample_s = pulseShaping_s_rx(1:D:end);
figure(6);
stem(downsample_s);
title('Downsampled Symbol for practical DAC');
