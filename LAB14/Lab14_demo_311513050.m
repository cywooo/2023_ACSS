clear all; clc; close all;

%% p1
signal_length = 1000;

for g = [1:10]
    nob_uni = g;
    DR_uni = 2^0; %+-
    uni = rand(1,signal_length)*2-1;
    code_map_uni = -DR_uni:((DR_uni*2)/(2^nob_uni)):(DR_uni-(DR_uni*2)/(2^nob_uni));
    uni_dec = F_point_decision(uni,code_map_uni);
    %noise_power(g) = 10*log10(mean(abs(uni-uni_dec).^2));
    SQNR_uni(g) = 10*log10(mean(uni.^2)/ mean(abs(uni-uni_dec).^2));
end
figure;
plot([1:length(SQNR_uni)],SQNR_uni);
xlabel("NOB");
title("DR = 1")

for g = [1:12]
    for k = [0:4]
        nob_nor = g;
        DR_nor = 2^k; %+-
        nor = randn(1,signal_length);
        code_map_nor = -DR_nor:((DR_nor*2)/(2^nob_nor)):(DR_nor-(DR_nor*2)/(2^nob_nor));
        nor_dec = F_point_decision(nor,code_map_nor);
        SQNR_nor(k+1,g) = 10*log10(norm(nor)/ norm(abs(nor-nor_dec)));
    end
end
figure;
plot([1:length(SQNR_nor(1,:))],SQNR_nor(1,:));
xlabel("NOB");
title("DR = 1")

figure;
plot([1:length(SQNR_nor(2,:))],SQNR_nor(2,:));
xlabel("NOB");
title("DR = 2")

figure;
plot([1:length(SQNR_nor(3,:))],SQNR_nor(3,:));
xlabel("NOB");
title("DR = 4")

figure;
plot([1:length(SQNR_nor(4,:))],SQNR_nor(4,:));
xlabel("NOB");
title("DR = 8")

%{
figure;
plot([1:length(uni)],uni);
hold on;
stem([1:length(uni_dec)],uni_dec,"O");
hold off;

figure;
plot([1:length(nor)],nor);
hold on;
stem([1:length(nor_dec)],nor_dec,"O");
hold off;
%}




%% p2

nob = 6;
DR = 2^3; %+-
code_map = -DR:((DR*2)/(2^nob)):(DR-(DR*2)/(2^nob));

signal_length = 2000;
s_in = randn(1,signal_length);

h = [1.5 1.8 0.8 0.6 ];
h = F_point_decision(h,code_map);

s1_noF = [s_in([1:end])*h(1)];
s2_noF = [0 s_in([1:end-1])*h(2)];
s3_noF = [0 0 s_in([1:end-2])*h(3)];
s4_noF = [0 0 0 s_in([1:end-3])*h(4)];
s_out_noF = s1_noF + s2_noF + s3_noF + s4_noF;

s1 = F_point_decision([s_in([1:end])*h(1)],code_map);
s2 = F_point_decision([0 s_in([1:end-1])*h(2)],code_map);
s3 = F_point_decision([0 0 s_in([1:end-2])*h(3)],code_map);
s4 = F_point_decision([0 0 0 s_in([1:end-3])*h(4)],code_map);

s5 =  F_point_decision((s1+s2),code_map);
s6 =  F_point_decision((s3+s5),code_map);
s_out_F_point =  F_point_decision((s4+s6),code_map);

figure; histogram(s_out_noF,2^nob);
title("p2");
figure; histogram(s_out_F_point,2^nob);
title("p3");

SQNR_P3 = 10*log10(mean(s_out_noF.^2)/ mean(abs(s_out_noF-s_out_F_point).^2));
