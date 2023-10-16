clear all
close all

%% Signal
N = 2^12;
x = randn(1,N);
h = SRRC_filter(2,2,0.3);

%% Direct Form
Sin = [zeros(1,length(h)-1) x];
% Floating point
for n=length(h):length(Sin)
    S1D(n) = h(1) * Sin(n);
    S2D(n) = h(2) * Sin(n-1);
    S3D(n) = S1D(n) + S2D(n);
    S4D(n) = h(3) * Sin(n-2);
    S5D(n) = S3D(n) + S4D(n);
    S6D(n) = h(4) * Sin(n-3);
    S7D(n) = S5D(n) + S6D(n);
    S8D(n) = h(5) * Sin(n-4);
    S9D(n) = S7D(n) + S8D(n);
    S10D(n) = h(6) * Sin(n-5);
    S11D(n) = S9D(n) + S10D(n);
    S12D(n) = h(7) * Sin(n-6);
    S13D(n) = S11D(n) + S12D(n);
    S14D(n) = h(8) * Sin(n-7);
    S15D(n) = S13D(n) + S14D(n);
    S16D(n) = h(9) * Sin(n-8);
    SoutD(n) = S15D(n) + S16D(n);
end

%plotFigure(Sin(9:end),S1D(9:end),S2D(9:end),S3D(9:end),S4D(9:end),S5D(9:end),S6D(9:end),S7D(9:end),S8D(9:end),S9D(9:end),S10D(9:end),S11D(9:end),S12D(9:end),S13D(9:end),S14D(9:end),S15D(9:end),S16D(9:end),SoutD(9:end));

% Fixed point
NOB_Dir = 9;

Sinq = zeros(1,length(Sin));
for n=1:length(Sin)
    Sinq(n) = ADC(Sin(n),NOB_Dir,4);
end

hq = zeros(1,length(h));
for n=1:length(h)
    hq(n) = ADC(h(n),NOB_Dir,1);
end

for n=length(hq):length(Sinq)
    S1Dq(n) = ADC(hq(1)*Sinq(n),NOB_Dir,0.2);
    S2Dq(n) = ADC(hq(2)*Sinq(n-1),NOB_Dir,0.4);
    S3Dq(n) = ADC(S1Dq(n)+S2Dq(n),NOB_Dir,0.4);
    S4Dq(n) = ADC(hq(3)*Sinq(n-2),NOB_Dir,0.2);
    S5Dq(n) = ADC(S3Dq(n)+S4Dq(n),NOB_Dir,0.5);
    S6Dq(n) = ADC(hq(4)*Sinq(n-3),NOB_Dir,1.5);
    S7Dq(n) = ADC(S5Dq(n)+S6Dq(n),NOB_Dir,1.5);
    S8Dq(n) = ADC(hq(5)*Sinq(n-4),NOB_Dir,2.5);
    S9Dq(n) = ADC(S7Dq(n)+S8Dq(n),NOB_Dir,3);
    S10Dq(n) = ADC(hq(6)*Sinq(n-5),NOB_Dir,1.5);
    S11Dq(n) = ADC(S9Dq(n)+S10Dq(n),NOB_Dir,4);
    S12Dq(n) = ADC(hq(7)*Sinq(n-6),NOB_Dir,0.2);
    S13Dq(n) = ADC(S11Dq(n)+S12Dq(n),NOB_Dir,4);
    S14Dq(n) = ADC(hq(8)*Sinq(n-7),NOB_Dir,0.4);
    S15Dq(n) = ADC(S13Dq(n)+S14Dq(n),NOB_Dir,4);
    S16Dq(n) = ADC(hq(9)*Sinq(n-8),NOB_Dir,0.2);
    SoutDq(n) = ADC(S15Dq(n)+S16Dq(n),NOB_Dir,4);
end

SQNR_D = 10*log10(mean(abs(SoutD(9:end)).^2)/mean(abs(SoutD(9:end)-SoutDq(9:end)).^2));
disp(['SQNR of Direct Form : ',num2str(SQNR_D),'dB']);

%% Transposed Form
Sin = [0 x];
% Floating point
S1T = zeros(1,length(Sin));
S3T = zeros(1,length(Sin));
S5T = zeros(1,length(Sin));
S7T = zeros(1,length(Sin));
S9T = zeros(1,length(Sin));
S11T = zeros(1,length(Sin));
S13T = zeros(1,length(Sin));
S15T = zeros(1,length(Sin));

for n=2:length(Sin)
    S1T(n) = h(9) * Sin(n);
    S2T(n) = h(8) * Sin(n);
    S3T(n) = S1T(n-1) + S2T(n);
    S4T(n) = h(7) * Sin(n);
    S5T(n) = S3T(n-1) + S4T(n);
    S6T(n) = h(6) * Sin(n);
    S7T(n) = S5T(n-1) + S6T(n);
    S8T(n) = h(5) * Sin(n);
    S9T(n) = S7T(n-1) + S8T(n);
    S10T(n) = h(4) * Sin(n);
    S11T(n) = S9T(n-1) + S10T(n);
    S12T(n) = h(3) * Sin(n);
    S13T(n) = S11T(n-1) + S12T(n);
    S14T(n) = h(2) * Sin(n);
    S15T(n) = S13T(n-1) + S14T(n);
    S16T(n) = h(1) * Sin(n);
    SoutT(n) = S15T(n-1) + S16T(n);
end
%plotFigure(Sin(2:end),S1T(2:end),S2T(2:end),S3T(2:end),S4T(2:end),S5T(2:end),S6T(2:end),S7T(2:end),S8T(2:end),S9T(2:end),S10T(2:end),S11T(2:end),S12T(2:end),S13T(2:end),S14T(2:end),S15T(2:end),S16T(2:end),SoutT(2:end));
    
% Fixed point
NOB_Tra = 9;

Sinq = zeros(1,length(Sin));
S1Tq = zeros(1,length(Sinq));
S3Tq = zeros(1,length(Sinq));
S5Tq = zeros(1,length(Sinq));
S7Tq = zeros(1,length(Sinq));
S9Tq = zeros(1,length(Sinq));
S11Tq = zeros(1,length(Sinq));
S13Tq = zeros(1,length(Sinq));
S15Tq = zeros(1,length(Sinq));
for n=1:length(Sin)
    Sinq(n) = ADC(Sin(n),NOB_Tra,4);
end

for n=2:length(Sinq)
    S1Tq(n) = ADC(hq(9)*Sinq(n),NOB_Tra,0.2);
    S2Tq(n) = ADC(hq(8)*Sinq(n),NOB_Tra,0.4);
    S3Tq(n) = ADC(S1Tq(n-1)+S2Tq(n),NOB_Tra,0.4);
    S4Tq(n) = ADC(hq(7)*Sinq(n),NOB_Tra,0.2);
    S5Tq(n) = ADC(S3Tq(n-1)+S4Tq(n),NOB_Tra,0.4);
    S6Tq(n) = ADC(hq(6)*Sinq(n),NOB_Tra,2);
    S7Tq(n) = ADC(S5Tq(n-1)+S6Tq(n),NOB_Tra,2);
    S8Tq(n) = ADC(hq(5)*Sinq(n),NOB_Tra,2.5);
    S9Tq(n) = ADC(S7Tq(n-1)+S8Tq(n),NOB_Tra,3);
    S10Tq(n) = ADC(hq(4)*Sinq(n),NOB_Tra,2);
    S11Tq(n) = ADC(S9Tq(n-1)+S10Tq(n),NOB_Tra,4);
    S12Tq(n) = ADC(hq(3)*Sinq(n),NOB_Tra,0.2);
    S13Tq(n) = ADC(S11Tq(n-1)+S12Tq(n),NOB_Tra,4);
    S14Tq(n) = ADC(hq(2)*Sinq(n),NOB_Tra,0.4);
    S15Tq(n) = ADC(S13Tq(n-1)+S14Tq(n),NOB_Tra,4);
    S16Tq(n) = ADC(hq(1)*Sinq(n),NOB_Tra,0.2);
    SoutTq(n) = ADC(S15Tq(n-1)+S16Tq(n),NOB_Tra,4);
end

SQNR_T = 10*log10(mean(abs(SoutT(2:end)).^2)/mean(abs(SoutT(2:end)-SoutTq(2:end)).^2));
disp(['SQNR of Transposed Form : ',num2str(SQNR_T),'dB']);

%% Hybrid Form
Sin = [zeros(1,4) x];
% Floating point
S9H = zeros(1,length(Sin));
S11H = zeros(1,length(Sin));
S13H = zeros(1,length(Sin));
S15H = zeros(1,length(Sin));

for n=5:length(Sin)
    S1H(n) = h(1) * Sin(n);
    S2H(n) = h(2) * Sin(n-1);
    S3H(n) = h(3) * Sin(n-1);
    S4H(n) = h(4) * Sin(n-2);
    S5H(n) = h(5) * Sin(n-2);
    S6H(n) = h(6) * Sin(n-3);
    S7H(n) = h(7) * Sin(n-3);
    S8H(n) = h(8) * Sin(n-4);
    S9H(n) = h(9) * Sin(n-4);
    S10H(n) = S8H(n) + S9H(n-1);
    S11H(n) = S7H(n) + S10H(n);
    S12H(n) = S6H(n) + S11H(n-1);
    S13H(n) = S5H(n) + S12H(n);
    S14H(n) = S4H(n) + S13H(n-1);
    S15H(n) = S3H(n) + S14H(n);
    S16H(n) = S2H(n) + S15H(n-1);
    SoutH(n) = S1H(n) + S16H(n);
end
%plotFigure(Sin(5:end),S1H(5:end),S2H(5:end),S3H(5:end),S4H(5:end),S5H(5:end),S6H(5:end),S7H(5:end),S8H(5:end),S9H(5:end),S10H(5:end),S11H(5:end),S12H(5:end),S13H(5:end),S14H(5:end),S15H(5:end),S16H(5:end),SoutH(5:end));

% Fixed point
NOB_Hyb = 9;

Sinq = zeros(1,length(Sin));
for n=1:length(Sinq)
    Sinq(n) = ADC(Sin(n),NOB_Hyb,4);
end
S9Hq = zeros(1,length(Sinq));
S11Hq = zeros(1,length(Sinq));
S13Hq = zeros(1,length(Sinq));
S15Hq = zeros(1,length(Sinq));
for n=5:length(Sinq)
    S1Hq(n) = ADC(hq(1)*Sinq(n),NOB_Hyb,0.2);
    S2Hq(n) = ADC(hq(2)*Sinq(n-1),NOB_Hyb,0.4);
    S3Hq(n) = ADC(hq(3)*Sinq(n-1),NOB_Hyb,0.2);
    S4Hq(n) = ADC(hq(4)*Sinq(n-2),NOB_Hyb,2);
    S5Hq(n) = ADC(hq(5)*Sinq(n-2),NOB_Hyb,3);
    S6Hq(n) = ADC(hq(6)*Sinq(n-3),NOB_Hyb,2);
    S7Hq(n) = ADC(hq(7)*Sinq(n-3),NOB_Hyb,0.2);
    S8Hq(n) = ADC(hq(8)*Sinq(n-4),NOB_Hyb,0.4);
    S9Hq(n) = ADC(hq(9)*Sinq(n-4),NOB_Hyb,0.2);
    S10Hq(n) = ADC(S8Hq(n)+S9Hq(n-1),NOB_Hyb,0.4);
    S11Hq(n) = ADC(S7Hq(n)+S10Hq(n),NOB_Hyb,0.4);
    S12Hq(n) = ADC(S6Hq(n)+S11Hq(n-1),NOB_Hyb,2);
    S13Hq(n) = ADC(S5Hq(n)+S12Hq(n),NOB_Hyb,3);
    S14Hq(n) = ADC(S4Hq(n)+S13Hq(n-1),NOB_Hyb,4);
    S15Hq(n) = ADC(S3Hq(n)+S14Hq(n),NOB_Hyb,4);
    S16Hq(n) = ADC(S2Hq(n)+S15Hq(n-1),NOB_Hyb,4);
    SoutHq(n) = ADC(S1Hq(n)+S16Hq(n),NOB_Hyb,4);
end

SQNR_H = 10*log10(mean(abs(SoutH(5:end)).^2)/mean(abs(SoutH(5:end)-SoutHq(5:end)).^2));
disp(['SQNR of Hybrid Form : ',num2str(SQNR_H),'dB']);

%% Plot Figure
figure
histogram(SoutD(9:end));
hold on
histogram(SoutT(2:end));
histogram(SoutH(5:end));
legend('Direct','Transposed','Hybrid');
title('Sout of three structures');

figure
subplot(3,1,1);
histogram(SoutD(9:end));
hold on
histogram(SoutDq(9:end));
legend('Original','Quanitized');
title('Direct Form');
subplot(3,1,2);
histogram(SoutT(2:end));
hold on
histogram(SoutTq(2:end));
legend('Original','Quanitized');
title('Transposed Form');
subplot(3,1,3);
histogram(SoutH(5:end));
hold on
histogram(SoutHq(5:end));
legend('Original','Quanitized');
title('Hybrid Form');

function p = plotFigure(Sin,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,Sout)
figure
subplot(3,6,1);
histogram(Sin);
title('Sin');
subplot(3,6,2);
histogram(S1);
title('S1');
subplot(3,6,3);
histogram(S2);
title('S2');
subplot(3,6,4);
histogram(S3);
title('S3');
subplot(3,6,5);
histogram(S4);
title('S4');
subplot(3,6,6);
histogram(S5);
title('S5');
subplot(3,6,7);
histogram(S6);
title('S6');
subplot(3,6,8);
histogram(S7);
title('S7');
subplot(3,6,9);
histogram(S8);
title('S8');
subplot(3,6,10);
histogram(S9);
title('S9');
subplot(3,6,11);
histogram(S10);
title('S10');
subplot(3,6,12);
histogram(S11);
title('S11');
subplot(3,6,13);
histogram(S12);
title('S12');
subplot(3,6,14);
histogram(S13);
title('S13');
subplot(3,6,15);
histogram(S14);
title('S14');
subplot(3,6,16);
histogram(S15);
title('S15');
subplot(3,6,17);
histogram(S16);
title('S16');
subplot(3,6,18);
histogram(Sout);
title('Sout');
end




    