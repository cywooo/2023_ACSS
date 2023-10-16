clear;close all;clc;
load('P3');

% -----------(b)-----------
uf = 64;
a = zeros(1,length(g3b_bit)*uf);
a(1:uf:end) = g3b_bit;
filt = SRRC(0.25,uf,10*uf);
r3b_x = conv(a,filt);
figure();plot(filt);
figure();plot(r3b_x);

% -----------(c)-----------
y3c = conv(r3b_x,filt);
r3c_y = y3c(1:uf:length(r3b_x));
figure();plot(y3c);
figure();stem(g3b_bit);hold on;stem(0.015*r3c_y(11:end));


% -----------(d)-----------
a = zeros(1,length(g3b_bit)*uf);
for i=1:length(g3b_bit)
    a(1+(i-1)*uf:i*uf) = g3b_bit(i);
end
r3d_x = conv(a,filt);
y3d = conv(r3b_x,filt);
snr_star = 0;
ph_star = 0;
for ph = 1:uf-1
    r3d_y = y3d(1+ph:uf:end);
    Esig = mean(abs(g3b_bit).^2);
    En = mean(abs(r3d_y(11:end)).^2)*0.0003-Esig;
    tmp = 10*log10(Esig/abs(En));
    if tmp>snr_star
        snr_star = tmp;
        ph_star = ph;
    end
end
phase = ph_star;
r3d_y = y3d(1+phase:uf:end);
Esig = mean(abs(g3b_bit).^2);
En = mean(abs(r3d_y(11:end)).^2)*0.0003-Esig;
snr_d = 10*log10(Esig/abs(En));
figure();stem(a);
figure();stem(g3b_bit);hold on;stem(0.015*r3d_y(11:end-10));
    
% -----------(e)-----------
a_e = zeros(1,length(g3b_bit)*uf);
a_e(1:uf:end) = g3b_bit;
f = linspace(-1/2,1/2,length(a_e));
figure();stem(f,abs(fftshift(fft(a_e))));
y_e = filter(r3e_IIRps,a_e);
yr_e = filter(r3e_IIRps,y_e);
r3e_y = yr_e(1:uf:end);
delay = 31;
En_e = mean(abs(r3e_y(1+delay:end)).^2)-Esig*0.01;
snr_e = 10*log10(Esig/abs(En_e));
figure();stem(g3b_bit);hold on;stem(100*r3e_y(1+delay:end));

% -----------(a)-----------
function srrc = SRRC(alpha,uf,span)
    T = uf;
    srrc = zeros(1,length(span));
    mp = 1/2*span;
    for i=1:1/2*span
        t=i+0.0001;
        srrc(mp+i) = (4*alpha/pi) * ( cos((1+alpha)*pi*t/T)+T*sin((1-alpha)*pi*t/T)/(4*alpha*t) ) / (1-(4*alpha*t/T)^2);
    end
    srrc(mp:-1:1) = srrc(mp+1:end);
end

