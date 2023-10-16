clear all
close all

N = 2^12;
x = randn(1,N);
Sin = [zeros(1,3) x];
%h = randn(1,4);
h = [-2.00793029130224 0.526057845795559 -0.741515913144450 -0.0861755628843450];

%% Practice 2
for n=4:length(Sin)
    S1(n) = h(1) * Sin(n);
    S2(n) = h(2) * Sin(n-1);
    S5(n) = S1(n) + S2(n);
    S3(n) = h(3) * Sin(n-2);
    S6(n) = S3(n) + S5(n);
    S4(n) = h(4) * Sin(n-3);
    Sout(n) = S4(n) + S6(n);
end

figure(1);
subplot(2,4,1);
histogram(Sin(4:end));
title('Sin');
subplot(2,4,2);
histogram(S1(4:end));
title('S1');
subplot(2,4,3);
histogram(S2(4:end));
title('S2');
subplot(2,4,4);
histogram(S3(4:end));
title('S3');
subplot(2,4,5);
histogram(S4(4:end));
title('S4');
subplot(2,4,6);
histogram(S5(4:end));
title('S5');
subplot(2,4,7);
histogram(S6(4:end));
title('S6');
subplot(2,4,8);
histogram(Sout(4:end));
title('Sout');

%% Practice3
hq = zeros(1,length(h));
for n=1:length(h)
    hq(n) = ADC(h(n),4,2);
end

Sin_q = zeros(1,length(Sin));
for n=1:length(Sin)
    Sin_q(n) = ADC(Sin(n),5,4);
end

for n=4:length(Sin_q)
    S1q(n) = ADC(hq(1)*Sin_q(n),8,8);
    S2q(n) = ADC(hq(2)*Sin_q(n-1),5,2);
    S5q(n) = ADC(S1q(n)+S2q(n),8,8);
    S3q(n) = ADC(hq(3)*Sin_q(n-2),7,4);
    S6q(n) = ADC(S3q(n)+S5q(n),8,8);
    S4q(n) = ADC(hq(4)*Sin_q(n-3),5,1);
    Soutq(n) = ADC(S4q(n)+S6q(n),8,8);
end

figure(2);
histogram(Sout);
hold on
histogram(Soutq);
title('Compare of Sout');

%% SQNR
SQNR = 10*log10(mean(abs(Sout).^2)/mean(abs(Sout-Soutq).^2))


