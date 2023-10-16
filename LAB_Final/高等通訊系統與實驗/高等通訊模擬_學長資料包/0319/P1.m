t=[1:128];
fs = 1/64;
fc = 1/4;

%windowsize = 6;
%b = (1/windowsize)*ones(1,windowsize);
b = poly([-0.9,-0.9]);
s1 = cos(2*pi*fs*t);
s2 = triangularPulse(1,32,63,t);

carrier = sqrt(2)*cos(2*pi*fc*t);

mod1=s1.*carrier;   %mod
dem1=mod1.*carrier;  %demo
mod2=s2.*carrier;   %mod
dem2=mod2.*carrier;  %demo

figure(1)
zplane(b,1);
figure(2)
freqz(b,1);
figure(3)
recover_s1=filter(b,1,dem1);
plot([s1' recover_s1'])
figure(4)
recover_s2=filter(b,1,dem2);
plot([s2' recover_s2'])

