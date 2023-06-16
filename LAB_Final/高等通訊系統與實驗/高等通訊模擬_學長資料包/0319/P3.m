t=[1:128];
fs = 1/64;
fc = 1/4;

% windowsize = 6;
% b = (1/windowsize)*ones(1,windowsize);
b = poly([-0.9,-0.9]);
           
rect = ones(1,32);

s1 = [zeros(1,32),ones(1,32),zeros(1,64)];
s2 = [zeros(1,32),conv(ones(1,32),ones(1,32)),zeros(1,33)];

e1 = sqrt(2)*exp(2*pi*fc*t);
e2 = sqrt(2)*exp(-2*pi*fc*t);

mod =(s1 + 1i*s2).*e1;
dem = mod.*e2;

recover=filter(b,1,dem);

figure(1)
zplane(b,1);
figure(2)
freqz(b,1);
figure(3)
r1 = real(recover);
plot([s1' r1'])

figure(4)
r2 = imag(recover);
plot([s2' r2'])

