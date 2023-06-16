clear all; close all; clc;
data = load('data_given');
a = data.a_4;

%% (a)
h = [1 -1.2 -0.25 0.3];
k = 20;

% syms z
% factor(1 - 1.2*z^(-1) - 0.25*z^(-2) + 0.3*z^(-3) )

eq1 = -(1.2).^(-10:-1); 
eq2 = (-0.5).^(0:5-1); 
eq3 = (0.5).^(0:5-1); 

eq_temp = conv(eq1,eq2);
eq = conv(eq_temp,eq3);
imp_4a = conv(eq,1);

h_imp = conv(h,1);
eqc_4a = conv(h_imp,eq);
delay = 10;
eqc_4a = eqc_4a(delay+1:end)
figure(1);
stem(imp_4a);
title('imp 4a')
figure(2);
stem(eqc_4a);
title('eqc 4a')

%% (b)
b = conv(h,a);

a_4hb = conv(eq,b);
delay = 10;
a_4hb = a_4hb(delay+1:end-10); %%%%

figure(3);
stem(a);
hold on 
stem(a_4hb);

e = a_4hb-a;
a_power = mean(abs(a).^2);
e_power = mean(abs(e).^2);
SNR_4b = 10*log10(a_power/e_power);

%% (c)
b = conv(h,a);

SNR = 10;
a_power = mean(abs(a).^2);
n_power = a_power/(10^(SNR/10));
n = sqrt(n_power)*randn(1,length(b));
sd_4c = std(n);
b_n = b+n;

a_4hc = conv(eq,b_n);
delay = 10;
a_4hc = a_4hc(delay+1:end-10); %%%%

figure(4)
stem(a)
hold on 
stem(a_4hc)

e = a_4hc-a;
a_power = mean(abs(a).^2);
e_power = mean(abs(e).^2);
SNR_4c = 10*log10(a_power/e_power);

%% (d)
detect = (a_4hc > 0)-(a_4hc < 0);
err = sum(detect ~= a);
SER_4d = err/length(a);

%% (e) 
b = conv(h,a);
SNR = 10;
a_power = mean(abs(a).^2);
n_power = a_power/(10^(SNR/10));
n = sqrt(n_power)*randn(1,length(b));
sd_4c = std(n);
b_n = b+n;

aeq = [-1:-1:-20];
a_4he = conv(aeq,b_n);
delay = 0;
a_4he = a_4hc(delay+1:end); 

figure(5)
stem(a)
hold on 
stem(a_4he)
legend('original','equlaizer')

detect = (a_4he > 0) - (a_4he < 0);
error = sum(detect ~= a);
SER_4e = error/length(a)