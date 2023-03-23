clear all; clc; close all;
%% LAB04 demo
n = [1:201];
fc = 1/4;
delay = 4;

w = linspace(-1/2,1/2,length(n));
x_I = 0 + 1 .*(n >= 50).*(n < 151);
x_Q = 0 + (n-50)/50 .*(n >= 50).*(n < 101)+(151-n)/50.*(n >= 101).*(n < 151);

% modulation
x = x_I.*(sqrt(2).* cos(2*pi*fc*n))+x_Q.*(-sqrt(2).* sin(2*pi*fc*n));
% channel effect
y = [zeros(1,delay) x([1:length(x)-delay])];
% demodulation
y_1 = y .* (sqrt(2).* cos(2*pi*fc*n));
y_2 = y .* (-sqrt(2).* sin(2*pi*fc*n));

% low-pass filter
h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=1.6*h(max(j-1,1)).*(j>1)-0.7300*h(max(j-2,1)).*(j>2)+...
        1.*(j == 0+1)+0.9.*(j == 1+1);       
end

[conv_y_1,conv_length] = conv_(y_1,h);
[conv_y_2,conv_length] = conv_(y_2,h);

figure(1);
subplot(4,2,1);
plot(n,x_I);
title("I-phase signal");
subplot(4,2,2);
plot(n,x_Q);
title("Q-phase signal");
subplot(4,2,3);
plot(n,y_1);
title("modulated y1");
subplot(4,2,4);
plot(n,y_2);
title("modulated y2");
subplot(4,2,5);
plot([1:conv_length],conv_y_1);
title("I-phase recovered by method 1");
subplot(4,2,6);
plot([1:conv_length],conv_y_2);
title("Q-phase recovered by method 1");


%%
m = x_I + x_Q*1i ;
carrier = exp(1i*2*pi*fc*n);
% modulation
x_2 = real( m *sqrt(2).* carrier);
% channel effect
y_2 = [zeros(1,delay) x_2([1:length(x_2)-delay])];
% demodulation
m_de = y_2 *sqrt(2).* (carrier.^-1);

%LPF
h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=1.6*h(max(j-1,1)).*(j>1)-0.7300*h(max(j-2,1)).*(j>2)+...
        1.*(j == 0+1)+0.9.*(j == 1+1);       
end
[conv_m_de,conv_length] = conv_(m_de,h);

x_I_de = real(conv_m_de);
x_Q_de = imag(conv_m_de);

figure(1);
subplot(4,2,7);
plot([1:conv_length],x_I_de);
title("I-phase recovered by method 2");
subplot(4,2,8);
plot([1:conv_length],x_Q_de);
title("Q-phase recovered by method 2");

