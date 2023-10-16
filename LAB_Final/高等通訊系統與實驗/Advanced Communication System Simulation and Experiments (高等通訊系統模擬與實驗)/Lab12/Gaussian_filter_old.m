function h = Gaussian_filter_old(M,BT,trun)
    n = [-M*trun:M*trun];
    h = exp((-2*(pi^2)/log(2))*((BT/M)^2).*(n.^2));
    C = sqrt(2*pi/log(2)) * BT;
    h = h*C;
    h = h/max(h);
end