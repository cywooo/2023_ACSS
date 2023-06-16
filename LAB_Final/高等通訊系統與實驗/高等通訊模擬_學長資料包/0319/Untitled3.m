n = 512;
fs = 2;
sd=0;
s=[1:n];
fq_1 = 1/64;
fq_2 = 1/4;
lpf=ones(1,fs)/fs;
y = sin(2*pi*fq_1*s);
x= triangularPulse(1,32,63,s)+ triangularPulse(65,96,127,s) + triangularPulse(129,160,191,s);
ca_1 = sqrt(2)*cos(2*pi*fq_2*s);
ca_2 = sqrt(2)*sin(2*pi*fq_2*s);
w=sd*randn(1,n);
mx=x.*ca_1;
mxn=mx+w;
cmx=mxn.*ca_1;
dmx=filter(lpf,1,cmx);
my=y.*ca_2;
myn=my+w;
cmy=myn.*ca_2;
mz = mx+my;
cmz = cmy + cmx;

dmzx=filter(lpf,0.995,cmx);
dmzy=filter(lpf,0.98,cmy);
figure(1)
 plot([dmzx' x'])
figure(2)
 plot([dmzy' y'])

