function up_s = up_sample(up,s)
    up_s = zeros(up,length(s)) ;
    up_s(1,:) = s;
    up_s = reshape(up_s,1,up*length(up_s));
end