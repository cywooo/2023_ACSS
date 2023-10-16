function SRRC = SRRC_filter(trun,M,alpha)
    n = [-M*trun:M*trun]+0.00001;
    nom = (4*alpha/pi) * (cos((1+alpha)*pi*n/M) + M*sin((1-alpha)*pi*n/M)/4/alpha./n);
    %nom = (4*alpha/pi) * (cos((1+alpha)*pi*n/M) + sinc((1-alpha)*n/M)/4/alpha*(1-alpha)*pi);
    den = 1-(4*alpha*n/M).^2;
    SRRC = nom./den;
end