clear all; close all; clc;
% signal
a = zeros(1,16);
a(1) = 1;

% channel
h = [1 -2 0.6];
b = conv(a,h);

% equalizer
const2 = 5/(2*sqrt(10)) * (-5+sqrt(10));
const1 = 5-const2;
w = [ const1 * -(1+(1/5)*sqrt(10)).^[-3:-1] const2 * (1-1/5*sqrt(10)).^[0:3] ];
a_hat = conv(b,w);
% normalize 
a_hat = a_hat/max(a_hat);

figure();stem(a_hat);

% plot figure
% figure();
% subplot(3,1,1);stem(a);hold on;
% subplot(3,1,2);stem(b);
% subplot(3,1,3);stem(a_hat);hold off;

% impulse response of equalizer
P1a_eq = conv(1,w);
figure();
stem(P1a_eq);

% equalized channel
P1a_imp = conv(conv(1,h),w);
figure();
stem(P1a_imp);

% SIR
Psig = 1;
Pif = (sum(abs(a_hat).^2) - 1)/length(a_hat);
P1b_SIR = 10*log10(Psig/Pif);