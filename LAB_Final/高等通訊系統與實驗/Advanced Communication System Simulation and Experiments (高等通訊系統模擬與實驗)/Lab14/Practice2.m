clear all
close all

N = 2^12;
x = randn(1,N);
Sin = [zeros(1,3) x];
h = randn(1,4);

for n=4:length(Sin)
    S1(n) = h(1) * Sin(n);
    S2(n) = h(2) * Sin(n-1);
    S5(n) = S1(n) + S2(n);
    S3(n) = h(3) * Sin(n-2);
    S6(n) = S3(n) + S5(n);
    S4(n) = h(4) * Sin(n-3);
    Sout(n) = S4(n) + S6(n);
end

figure(1);
subplot(2,4,1);
histogram(Sin(4:end));
title('Sin');
subplot(2,4,2);
histogram(S1(4:end));
title('S1');
subplot(2,4,3);
histogram(S2(4:end));
title('S2');
subplot(2,4,4);
histogram(S3(4:end));
title('S3');
subplot(2,4,5);
histogram(S4(4:end));
title('S4');
subplot(2,4,6);
histogram(S5(4:end));
title('S5');
subplot(2,4,7);
histogram(S6(4:end));
title('S6');
subplot(2,4,8);
histogram(Sout(4:end));
title('Sout');

% figure(2);
% histogram(h);
% title('h');
