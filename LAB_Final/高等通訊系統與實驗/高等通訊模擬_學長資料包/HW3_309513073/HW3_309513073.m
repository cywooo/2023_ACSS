
%_______________________FIR______________________________________
% a = poly([0.9+0.41j,0.59+0.8j,0.9-0.41j,0.59-0.8j]);       %zero numerator���l
% b = 1;%pole 
%_______________________IIR__________________
a = poly([0.9+0.41j,0.59+0.8j,0.9-0.41j,0.59-0.8j]);                        %zero numerator���l
b = poly([-0.07+0.8j, -0.2+0.323j, -0.2-0.323j, -0.07-0.8j]);         %pole ���� 
figure(1)
zplane(a,b);
figure(2)
freqz(a,b,100);
t = 1:128;
s1 = cos(2*pi*0.25*t);
s2 = 0.2*cos(2*pi*0.125*t);   %110 17 �C�W  =���Tv
s = s1+s2;

figure(3)     %���T�e���T��
sf=abs(fft(s));
stem(sf)


figure(4)       % �T���M���T�@�_�o�i���T��
S = filter(a,b,s); % ���l(zeros) , ����(pole)
Sf=abs(fft(S));
stem(Sf);   

SNR = 10*log10(mean(abs(S).^2)/mean(abs(s2).^2));

