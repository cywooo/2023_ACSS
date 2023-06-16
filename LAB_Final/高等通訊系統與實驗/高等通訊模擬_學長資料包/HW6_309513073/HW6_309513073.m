F=7;

df=F;
uf=F;


ad = a1(1:df:end); % down

Ad = [ad;zeros(size(a1,1)-size(ad,1),1)];
n = size(ad,1);

adu = zeros(n*uf,1);
adu(1:uf:end)=ad;

fadu = filter(Hd,adu); %filter


figure(1)
stem(fftshift(abs(fft(a1))))
title('FFT of signal')
figure(2)%¿W√–≈‹ºe downsampling
stem(fftshift(abs(fft(Ad))))
title('downsampling')
figure(3)
stem(fftshift(abs(fft(adu))))
title('upsampling')
figure(4)% upsampling
stem(1.3.*fftshift(abs(fft(fadu))))
title('filter')
figure(5)
subplot(2,1,1),plot(a1)
title('origin signal')
subplot(2,1,2),plot(adu(35:end))
title('recover signal')
