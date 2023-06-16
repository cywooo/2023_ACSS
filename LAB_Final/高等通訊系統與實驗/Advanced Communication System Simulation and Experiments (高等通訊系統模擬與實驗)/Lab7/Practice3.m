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

%% Design an IIR DMA filter for a upsampling process with L=32
% Upsampling
L = 32;
upsample_sym = zeros(1,N*L);
upsample_sym(1:L:end) = symbol;

% Pulse shaping (IIR DMA)
pulseShaping_s = filter(IIR_DMA_filter,upsample_sym); % Call IIR_DMA_filter function
figure(2);
stem(pulseShaping_s);
title('Pulse Shaping Symbol for IIR DMA');

%% Let the same filter be used at the receiver
pulseShaping_s_rx = L*filter(IIR_DMA_filter,pulseShaping_s);
figure(3);
stem(pulseShaping_s_rx);
title('Received Pulse Shaping Symbol for IIR DMA');

%% Conduct detection at the receiver
D = 32;
pulseShaping_s_rx = pulseShaping_s_rx(51:end); %25*2
downsample_s = pulseShaping_s_rx(1:D:end);
figure(4);
stem(downsample_s);
title('Downsampled Symbol for IIR DMA');
