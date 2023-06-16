sn = 50;
b = randi([0,1],1,sn);
b(b==0)=-1;

%SRRC
% a = 0.25;
% M = 4;
% n = [-128:127]-0.001;
% SRRC = (4*a/pi) *( cos((1+a)*pi*n./M) + (M*sin((1-a)*pi*n./M) ./(4*a*n))  )./ (1-(4*a*n./M).^2 ) ;

h = rcosdesign(0.25,10,4,'sqrt');
%up sampling 4
uf1 = 4;
bup = zeros(1,uf1*sn);
bup(1:uf1:end) = b;
fbup = conv(bup,h,'same');

%up sampling total = 16
uf2 = 4;
fbupup = zeros(1,uf2*size(fbup,2));
fbupup(1:uf2:end) = fbup;
ffbupup = conv(fbupup,h,'same');

figure(1)
subplot(5,1,1);stem(abs(fft(b)));
subplot(5,1,2);stem(abs(fft(bup)));
subplot(5,1,3);stem(abs(fft(fbup)));
subplot(5,1,4);stem(abs(fft(fbupup)));
subplot(5,1,5);stem(abs(fft(ffbupup)));
%% p2
np = 0;
r = ffbupup + np*randn(1,length(ffbupup));
%Receiver
rf = conv(r,h,'same');
rfd = rf(1:uf2:end);
rfdf = conv(rfd,h,'same');
rfdfd = rfdf(1:uf1:end);

figure(6)
stem([b(1:end)' rfdfd(1:end)']);
%% p3
IIRfbupup = filter(IIR,fbupup);
r2 = IIRfbupup + np*randn(1,length(IIRfbupup));
rf2 = filter(IIR,r2);
rfd2 = rf2(1:uf2:end);
rfdf2 = conv(rfd2,h,'same');
rfdfd2 = rfd2(1:uf1:end);

figure(7)
stem(abs(fft(IIRfbupup)));
figure(8)
stem([b(1:end-1)' 8*rfdfd2(2:end)']);
