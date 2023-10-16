function RC = RC_filter(trun,M,alpha)
    n = [-M*trun:M*trun]+0.00001;
    RC = sinc(n/M) .* cos(2*pi*alpha*n/2/M) ./ (1 - 16*alpha^2*n.^2/4/M^2);
end