clear all;
n=40;
x(1:n)=1;
y=conv(x,x);
np=mean(y.^2)*10^(-1.5);
v = sqrt(np)*randn(1,n*2-1);
s = y + v;
SNR_=10*log10(mean(y.^2)/mean(v.^2));

for w = 1:n
    sf = fft(s);
    sfw = zeros(1,2*n-1);
    sfw(1:w)=sf(1:w);
    sfw(2*n+1-w:2*n-1)=sf(2*n+1-w:2*n-1);
    r = ifft(sfw);
    MSE(w) = mean((r-y).^2);%immse(r,y);
end

[M,I]=min(MSE);
figure(1)
plot([r' y'])
figure(2)
stem(abs(sfw))


