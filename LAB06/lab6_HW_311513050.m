clear all; clc; close all;
%% parameters
fs = 8000;
load("A1.mat");
a1 = a1'; % sound
n = [1:length(a1)]/fs; % timeline

s = a1;
S = fftshift(fft(s));
% sound(s,fs); listen to the data

figure(1);
subplot(3,1,1);
plot(n,s);
title("s");
subplot(3,1,[2:3]);
plot(linspace(-1/2,1/2,length(S)),abs(S));
title("|S|");

%% HW part01
% For a given signal,
% try to conduct downsamping with a largest factor without causing distortion.

% show different down sample rate
% to see what will happen 

down_test = [1:6];

for k = [1:length(down_test)]
    down_test_s = down_sample(down_test(k),s);
    down_test_S = fftshift(fft(down_test_s));
    
    figure(2);
    plot(linspace(-1/2,1/2,length(down_test_S)),abs(down_test_S));
    if k == 1
        hold on; 
    end
end
hold off;
legend("down = 1","down = 2","down = 3","down = 4","down = 5","down = 6",'FontSize',12);
%%  HW part02
% choose downsample rate 3
% upsample 3 for recovering
down = 3;
up = 3;

down_s = down_sample(down,s)*down;
down_S = fftshift(fft(down_s));
%sound(down_s,fs/down);

up_s = up_sample(up,down_s);
up_S = fftshift(fft(up_s));
%sound(up_s,fs/down*up);

LPF_s = filter(lab6_HW_LPF,up_s);
LPF_S = fftshift(fft(LPF_s));
%sound(LPF_s,fs/down*up);


figure(3);
plot([1:length(s)],s);
hold on;
plot([1:length(LPF_s)],LPF_s);
hold off;
title("Recover comparing",'FontSize',12);
legend("Original","Recovering",'FontSize',12)

figure(4);
plot(linspace(-1/2,1/2,length(LPF_S)),abs(LPF_S));
hold on;
plot(linspace(-1/2,1/2,length(S)),abs(S));
hold off;
title("Recover comparing",'FontSize',12);
legend("Recovering |S|","Original |S|",'FontSize',12)

figure(5);
subplot(4,1,1);
plot(linspace(-1/2,1/2,length(S)),abs(S));
title_type = "Oringinal |S|";
title(title_type,'FontSize',12);
subplot(4,1,2);
plot(linspace(-1/2,1/2,length(down_S)),abs(down_S));
title_type = "|S| downsample by "+ num2str(down);
title(title_type,'FontSize',12);
subplot(4,1,3);
plot(linspace(-1/2,1/2,length(up_S)),abs(up_S));
title_type = "|S_d_o_w_n| upsample by "+ num2str(up);
title(title_type,'FontSize',12);
subplot(4,1,4);
plot(linspace(-1/2,1/2,length(LPF_S)),abs(LPF_S));
title_type = "|S_u_p| pass through LPF";
title(title_type,'FontSize',12);
