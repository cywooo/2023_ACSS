close all;

M=16;
k=log2(M);
S = 60000;
x = randi([0 1],S,k) ;

map = zeros(1,S);
ymap = zeros(1,S);
Iaxis = zeros(1,S);
Qaxis = zeros(1,S);
I2axis = zeros(1,S);
Q2axis = zeros(1,S);
for i = 1:S
   if x(i,1:2)== [0 0]
       Iaxis(i)=-3;
   elseif x(i,1:2)== [0 1]
       Iaxis(i)=-1;
   elseif x(i,1:2)== [1 1]
       Iaxis(i)=1;
   elseif x(i,1:2)== [1 0]
       Iaxis(i)=3;
   end
   if x(i,3:4)== [0 0]
       Qaxis(i)=-3j;
   elseif x(i,3:4)== [0 1]
       Qaxis(i)=-1j;
   elseif x(i,3:4)== [1 1]
       Qaxis(i)=1j;
   elseif x(i,3:4)== [1 0]
       Qaxis(i)=3j;
   end
   map(i)=Iaxis(i)+Qaxis(i);
end





snr = [5:15];
SNR = 10.^(snr./10);
SER  = zeros(size(snr));
TSER  = zeros(size(snr));
err = zeros(size(snr));
t = 1;
loop = 1;
for k = snr
    Es = sum(abs(map).^2)/S;
    N0 = Es/10^(k/10);
   
    for j = 1:loop
        n = (1/sqrt(2))*(sqrt(N0)*randn(1,S) +sqrt(N0)*randn(1,S)*1j);
        y = map + n ;
        
        for i = 1:S
           if  real(y(i)) < -2
               I2axis(i)=-3;
           elseif -2< real(y(i)) && real(y(i)) < 0
               I2axis(i)=-1;
           elseif 0< real(y(i))  && real(y(i)) < 2
               I2axis(i)=1;
           elseif 2< real(y(i)) 
               I2axis(i)=3;
           end
           if  imag(y(i)) < -2
               Q2axis(i)=-3j;
           elseif -2< imag(y(i)) && imag(y(i)) < 0
               Q2axis(i)=-1j;
           elseif 0< imag(y(i)) && imag(y(i))< 2
               Q2axis(i)=1j;
           elseif 2< imag(y(i)) 
               Q2axis(i)=3j;
           end
           ymap(i)=I2axis(i)+Q2axis(i);

           test = ymap(i) - map(i) ; 
           if test ~=0
               err(t) = err(t) +1;
           end
        end

    end
    
    SER(t) = (err(t)/S)/loop;
    TSER(t) = 2*(3/4) * erfc(sqrt(SNR(t)/10))-((3/4) * erfc(sqrt(SNR(t)/10)))^2;
    t = t+1;
 
end

scatterplot(y);
figure()
semilogy(snr,SER, snr,TSER )
legend('Empirical SER','Theoretical SER');
xlabel('SNR(dB)');
ylabel('Symbol Error Rate');

