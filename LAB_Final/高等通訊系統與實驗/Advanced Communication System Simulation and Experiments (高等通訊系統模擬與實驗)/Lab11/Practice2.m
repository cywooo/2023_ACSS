clear all
close all

% Channel
Ldma = 32;
tao = Ldma-1;
h = [1 zeros(1,tao) -0.5];

% Inverse response
ZF_equalizer = 0.5.^(0:Ldma-1);

% Convolution
result = conv(h,ZF_equalizer);

figure(1);
stem(result);