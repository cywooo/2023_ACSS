t=[1:128];
fs = 1/64;
fc = 1/4;
windowsize = 20;
b = (1/windowsize)*ones(1,windowsize);
% b = poly([-0.9,-0.9]);
           
rect = ones(1,32);
s1 = [zeros(1,32),ones(1,32),zeros(1,64)];
s2 = [zeros(1,32),conv(ones(1,32),ones(1,32)),zeros(1,33)];


dt = 0.008;
IC = sqrt(2)*cos(2*pi*fc*(t-dt));
QC = -sqrt(2)*sin(2*pi*fc*(t-dt));

mod1 = s1.*IC;  
mod2 = s2.*QC;
mod = mod1 + mod2;
dem1 = mod.*IC;
dem2 = mod.*QC;


e1 = sqrt(2)*exp(2*pi*fc*(t-dt));
e2 = sqrt(2)*exp(-2*pi*fc*(t-dt));

modB =(s1 + 1i*s2).*e1;
demB = modB.*e2;

at = t-dt;

figure(1)
recover_s1=filter(b,1,dem1);
plot([s1' recover_s1'])
figure(2)
recover_s2=filter(b,1,dem2);
plot([s2' recover_s2'])


recoverB = filter(b,1,demB);
figure(3)
r1 = real(recoverB);
plot([s1' r1'])
figure(4)
r2 = imag(recoverB);
plot([s2' r2'])



