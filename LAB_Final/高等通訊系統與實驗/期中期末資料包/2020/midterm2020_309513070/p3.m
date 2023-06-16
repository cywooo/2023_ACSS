clear all; close all; clc;
load('G_data.mat');
figure();plot(x3);
figure();plot(h3);

x3_e = [zeros(1,length(h3)-1) x3 zeros(1,length(h3)-1)];
y3 = zeros(1,length(x3)+length(h3)-1);
for i = 1:length(y3)
    for j=1:length(h3)
        y3(i) = y3(i) + x3_e(i+j-1)*h3(end-j+1);
    end
end
figure();plot(y3);hold on;plot(conv(x3,h3));