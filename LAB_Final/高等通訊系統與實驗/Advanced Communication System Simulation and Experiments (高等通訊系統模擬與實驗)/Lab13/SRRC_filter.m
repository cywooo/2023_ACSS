function SRRC = SRRC_filter(trun,M,alpha)
    n = [-M*trun:M*trun]+0.00001;
    h_nom = (4*alpha/pi) * (cos((1+alpha)*pi*n/M) + M*sin((1-alpha)*pi*n/M)./(4*alpha*n));
    h_den = 1 - (4*alpha*n/M).^2;
    SRRC = h_nom ./ h_den;
    SRRC = SRRC./sqrt(sum(abs(SRRC).^2));
end