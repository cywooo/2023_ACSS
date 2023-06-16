function h = Gaussian_filter(BT,M,trun)
    n = [-M*trun:M*trun];
    C = sqrt(2*pi/log(2))*BT;
    h = C * exp(((-2*(pi)^2)/log(2)) * ((BT/M)^2) .* (n.^2));
    h = h ./ sqrt(sum(abs(h).^2));
end
