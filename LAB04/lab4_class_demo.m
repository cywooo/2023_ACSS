clear all; clc; close all;
%% LAB04 demo
n = [1:201];
fc = 1/4;

w = linspace(-1/2,1/2,length(n));
x_I = 0 + 1 .*(n >= 50).*(n < 151);
%x_I = repmat([0 0 0 0 0 1 1 1 0 0 0 0 0] , [1,51]);
%x_I = x_I(1:length(n));
x_Q = 0 + (n-50)/50 .*(n >= 50).*(n < 101)+(151-n)/50.*(n >= 101).*(n < 151);
%x_Q = repmat([0 0 0 0.25 0.5 0.75 1 0.75 0.5 0.25 0 0 0] , [1,51]);
%x_Q = x_Q(1:length(n));

carrier = exp(1i*2*pi*fc*n);

figure(1);
subplot(2,1,1);
plot(n,x_I);
title("I");
subplot(2,1,2);
plot(n,x_Q);
title("Q");

m = x_I + x_Q*1i ;


m_ch = real( m *sqrt(2).* carrier);
m_de = m_ch *sqrt(2).* (carrier.^-1);

%LPF
LPF = 0 + 2 .*(w > -0.45).*(w < 0.45);

h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=1.6*h(max(j-1,1)).*(j>1)-0.7300*h(max(j-2,1)).*(j>2)+...
        1.*(j == 0+1)+0.9.*(j == 1+1);       
end
[conv_m_de,conv_length] = conv_(m_de,h);

M = fftshift(fft(m));
M_ch = fftshift(fft(m_ch));
M_de = fftshift(fft(m_de));
M_LPF = M_de .* LPF;


figure(2);
subplot(5,1,1);
plot(w,abs(M));
title("|M|");
subplot(5,1,2);
plot(w,abs(M_ch));
title("|M ch|");
subplot(5,1,3);
plot(w,abs(M_de));
title("|M de|");
subplot(5,1,4);
plot(w,LPF);
title("LPF");
subplot(5,1,5);
plot(w,abs(M_LPF));
title("|M LPF|");

m_LPF = ifft(fftshift(M_LPF));

x_I_de = real(conv_m_de);
x_Q_de = imag(conv_m_de);

figure(3);
subplot(2,1,1);
plot([1:conv_length],abs(x_I_de));
title("I");
subplot(2,1,2);
plot([1:conv_length],abs(x_Q_de));
title("Q");


figure(4);
subplot(4,1,1);
plot(n,abs(m));
title("m");
subplot(4,1,2);
plot(n,m_ch);
title("m ch");
subplot(4,1,3);
plot(n,abs(m_de));
title("m de");
subplot(4,1,4);
plot(n,abs(m_LPF));
title("m LPF");

%% 

x = x_I.*(sqrt(2).* cos(2*pi*fc*n))+x_Q.*(-sqrt(2).* sin(2*pi*fc*n));
y = x;

y_1 = y .* (sqrt(2).* cos(2*pi*fc*n));
y_2 = y .* (-sqrt(2).* sin(2*pi*fc*n));

h=zeros(1,length(n));
for j = 1:length(n)
    h(j)=1.6*h(max(j-1,1)).*(j>1)-0.7300*h(max(j-2,1)).*(j>2)+...
        1.*(j == 0+1)+0.9.*(j == 1+1);       
end

[conv_y_1,conv_length] = conv_(y_1,h);
[conv_y_2,] = conv_(y_2,h);

figure(5);
subplot(4,1,1);
plot(n,y_1);
subplot(4,1,2);
plot(n,y_2);
subplot(4,1,3);
plot([1:conv_length],conv_y_1);
subplot(4,1,4);
plot([1:conv_length],conv_y_2);



