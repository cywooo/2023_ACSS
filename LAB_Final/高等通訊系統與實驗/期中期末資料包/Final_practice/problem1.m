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
xm1_1 = mmse*y


detect1 = [];
detect2 = [];

for i = 1:length(xm1_1(1,:))
    if (real(xm1_1(1,i)) > 0 && imag(xm1_1(1,i))> 0 )
        detect1(i) = 1+1j;
    elseif (real(xm1_1(1,i)) > 0 && imag(xm1_1(1,i)) < 0 )
        detect1(i) = 1-1j;
    elseif  ( real(xm1_1(1,i)) < 0 && imag(xm1_1(1,i)) > 0 ) 
        detect1(i) = -1+1j;
    elseif ( real(xm1_1(1,i)) < 0 && imag(xm1_1(1,i)) < 0 ) 
        detect1(i) = -1-1j;
    end
end

for i = 1:length(xm1_1(2,:))
    if (real(xm1_1(2,i)) > 0 && imag(xm1_1(2,i))> 0 )
        detect2(i) = 1+1j;
    elseif (real(xm1_1(2,i)) > 0 && imag(xm1_1(2,i)) < 0 )
        detect2(i) = 1-1j;
    elseif  ( real(xm1_1(2,i)) < 0 && imag(xm1_1(2,i)) > 0 ) 
        detect2(i) = -1+1j;
    elseif ( real(xm1_1(2,i)) < 0 && imag(xm1_1(2,i)) < 0 ) 
        detect2(i) = -1-1j;
    end
end

error1 = sum(x(1,:)~=detect1);
error2 = sum(x(2,:)~=detect2);
SER_1 = (error1+error2)/(length(x)*2)


