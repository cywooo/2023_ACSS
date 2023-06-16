t=[1:128];
fs = 1/64;
fc = 1/4;

% windowsize = 6;
% b = (1/windowsize)*ones(1,windowsize);
b = poly([-0.9,-0.9]);
           
rect = ones(1,32);

s1 = [zeros(1,32),ones(1,32),zeros(1,64)];
s2 = [zeros(1,32),conv(ones(1,32),ones(1,32)),zeros(1,33)];


IC = sqrt(2)*cos(2*pi*fc*t);
QC = -sqrt(2)*sin(2*pi*fc*(t));

mod1 = s1.*IC;
mod2 = s2.*QC;


mod = mod1 + mod2;
 
modd = circshift(mod,12);
dem1 = modd.*IC;
dem2 = modd.*QC;


figure(3)
recover_s1=filter(b,1,dem1);
plot([s1' recover_s1'])
figure(4)
recover_s2=filter(b,1,dem2);
plot([s2' recover_s2'])



