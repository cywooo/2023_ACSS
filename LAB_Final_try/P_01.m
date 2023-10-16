clear all; close all; clc;
data = load('data_given');
x = data.x_1;
y = data.y_1;
H = data.H_1;
w = y-H*x;

x1_power = mean(abs(x(1,:)).^2);
n1_power = mean(abs(w(1,:)).^2);

x2_power = mean(abs(x(2,:)).^2);
n2_power = mean(abs(w(2,:)).^2);

% ZF detector
ZF = pinv(H)*y;

% MMSE
rho = [x1_power/n1_power  x2_power/n2_power];
mmse = pinv(H'*H + diag(pinv(rho)))*H'; %%%%%
xm1_1 = mmse*y;
xm1_1_detect = sign(real(xm1_1)) + sign(imag(xm1_1))*1i;

SER_1 = sum(xm1_1_detect~= x,"all")/(length(x)*2);