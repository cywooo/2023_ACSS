
%_______________________IIR______________________________________
a = poly([-0.1+0.6j,-0.1-0.6j,-0.15+0.56j,-0.15-0.56j]);                        %zero numerator分子
b = poly([0.35+0.5i, 0.35-0.5j, 0.52+0.47j, 0.52-0.47j]);%pole 分母
figure(1)
zplane(a,b);
figure(2)
freqz(a,b,100);
t = 1:125;
s1 = cos(2*pi*0.25*t);
s2 = cos(2*pi*0.125*t);   %110 17 低頻
s = s1+s2;
v = 0.05*randn(1,125);
y = s + v;


figure(3)     %雜訊前的訊號
sf=abs(fft(s));
stem(sf)

figure(4)     %加上雜訊後的訊號
yf=abs(fft(y));
stem(yf)

figure(5)       % 訊號和雜訊一起濾波的訊號
Y = filter(b,a,y);
Yf=fft(Y);
stem(abs(Yf));   

H=sum(a)/(1+sum(b));
G = 1/H;   % for normalized Gain 的 para
e = y-s;
SNR = 10*log10(mean(abs(s*G).^2)/mean(abs(e).^2));

%_______________________FIR______________________________________
% a = poly([-0.1+0.6j,-0.1-0.6j,-0.15+0.56j,-0.15-0.56j]);    %zero numerator分子
% b = 1;%pole 
% figure(1)
% zplane(a,b);
% figure(2)
% freqz(a,b,100);
% t = 1:125;
% s1 = cos(2*pi*0.25*t);
% s2 = cos(2*pi*0.125*t);   %110 17 低頻
% s = s1+s2;
% v = 0.08*randn(1,125);
% y = s + v;
% 
% 
% figure(3)     %雜訊前的訊號
% sf=abs(fft(s));
% stem(sf)
% 
% figure(4)     %加上雜訊後的訊號
% yf=abs(fft(y));
% stem(yf)
% 
% figure(5)       % 訊號和雜訊一起濾波的訊號
% Y = filter(b,a,y);
% Yf=fft(Y);
% stem(abs(Yf));   
% 
% H=sum(a)/(1+sum(b));
% G = 1/H;   % for normalized Gain 的 para
% e = y-s;
% SNR = 10*log10(mean(abs(s*G).^2)/mean(abs(e).^2));

