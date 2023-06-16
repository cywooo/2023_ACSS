
%_______________________FIR______________________________________
% a = poly([0.9+0.41j,0.59+0.8j,0.9-0.41j,0.59-0.8j]);       %zero numerator分子
% b = 1;%pole 
%_______________________IIR__________________
a = poly([0.9+0.41j,0.59+0.8j,0.9-0.41j,0.59-0.8j]);                        %zero numerator分子
b = poly([-0.07+0.8j, -0.2+0.323j, -0.2-0.323j, -0.07-0.8j]);         %pole 分母 
figure(1)
zplane(a,b);
figure(2)
freqz(a,b,100);
t = 1:128;
s1 = cos(2*pi*0.25*t);
s2 = 0.2*cos(2*pi*0.125*t);   %110 17 低頻  =雜訊v
s = s1+s2;

figure(3)     %雜訊前的訊號
sf=abs(fft(s));
stem(sf)


figure(4)       % 訊號和雜訊一起濾波的訊號
S = filter(a,b,s); % 分子(zeros) , 分母(pole)
Sf=abs(fft(S));
stem(Sf);   

SNR = 10*log10(mean(abs(S).^2)/mean(abs(s2).^2));

