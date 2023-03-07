function [c,out_length]  = conv_(a,b)
 out_length = length(a)+length(b)-1;
    c = zeros(1,out_length);
    min_length = length(a)*(length(a)<length(b)) + length(b)*(length(b)<=length(a));
    for i = 1:out_length
        overlap_length = min([i,min_length,-1*(i-out_length-1)]);
        a_index = -(-1*max([-i,-length(a)])-overlap_length+1):-1:max([-i,-length(a)]);
        b_index = max([-i,-length(b)]):1:-(-1*max([-i,-length(b)])-overlap_length+1);
        c(out_length-i+1) = sum(a(a_index+length(a)+1).* b(b_index+length(b)+1));
    end 
end