
load('C:\Users\USER\AppData\Local\Temp\Rar$DRa21248.35195\G_data.mat')

%%
y1 = abs(fft(x1));
y2 = fftshift(abs(fft(x1)));
x = (-1/2:1/2) /4096;
xlabel=(-4096/2:4096/2-1)*(10^6/4096);
figure(1)
stem(xlabel,y2);
freq1 = (find(y1>0.3) / 4096) *10^6  ;
freq2 = (xlabel(find(y2>0.3) ))  ;

%%

a = [zeros(1,64) ones(1,1) zeros(1,64)];
af = filter(xf2,1,a);

p = phase(xf2);


figure(2)
% stem(fftshift(abs(xf2)));
plot(p);
% figure(2)
figure(7)
stem([15*a' abs(af')]);

%%

y3 = zeros(1,length(x3)+length(h3)-1);
Y3 = conv(x3,h3);        
YY3 = filter(h3,1,x3);   
for i = 1:length(x3)   
   y3 = y3 + [zeros(1,i-1), x3(i)*h3, zeros(1,length(y3)-length(h3)-i+1)]; 
    
end

%%
load filter4b.mat
load filter4c.mat
n = length(x4);
SPW = (mean(x4.^2) - mean(x4)^2 );
SNR = 10;
NPW = SPW * 10^(-SNR/10);    %noise power

xn4a = sqrt(NPW)*randn(1,n);    
sd4a = x4 + xn4a;
SNR =10*log10(mean(sd4a.^2 - mean(sd4a)^2)/mean(xn4a.^2));

y4b = filter(filter4b,sd4a);   % y
xnf4b = filter(filter4b,x4);        % sb
e4b = y4b - xnf4b;
snr4b =10*log10(mean(xnf4b.^2 - mean(xnf4b)^2)/mean(e4b.^2));
figure(4)
stem(abs(fft(xnf4b)));

y4c = filter(filter4c,sd4a);   % y
xnf4c = filter(filter4c,x4);        % sb
e4c = y4c - xnf4c;
snr4c =10*log10(mean(xnf4c.^2 - mean(xnf4c)^2)/mean(e4c.^2));


%%
% M = uf;
% n = [-128:127]-0.0001;
% SRRC = (4*a/pi)* ( cos((1+a)*pi*n./M ) + M*sin((1-a)*pi*n./M)  ./ (4*a*n)  ) ./ (1-(4*a*n./M).^2);




















































