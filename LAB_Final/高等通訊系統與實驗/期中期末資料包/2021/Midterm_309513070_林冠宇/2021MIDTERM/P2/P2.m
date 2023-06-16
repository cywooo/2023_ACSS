clear;close all;clc;
load('P2');

r2_mod = zeros(1,length(g2_x)/6);
k=1;
for i=1:6:length(g2_x)
    bit_seq = g2_x(i:i+5);
    for j=1:2
        idx = 1+3*(j-1);
        if bit_seq(idx:idx+2) == [0 0 0]
            val = 7;
        elseif bit_seq(idx:idx+2) == [0 0 1]
            val = 5;
        elseif bit_seq(idx:idx+2) == [0 1 1]
            val = 3;
        elseif bit_seq(idx:idx+2) == [0 1 0]
            val = 1;
        elseif bit_seq(idx:idx+2) == [1 1 0]
            val = -1;
        elseif bit_seq(idx:idx+2) == [1 1 1]
            val = -3;
        elseif bit_seq(idx:idx+2) == [1 0 1]
            val = -5;
        else
            val = -7;
        end
        if j==1
            r2_mod(k) = r2_mod(k) + val;
        else
            r2_mod(k) = r2_mod(k) + 1j*val;
        end
    end
    k=k+1;
end