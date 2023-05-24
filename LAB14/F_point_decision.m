function s_out = F_point_decision(s_in,map)
    divided = abs(ones(length(map),1)*s_in - map');
    min_div = min(divided);
    [s_x s_y] = size(s_in);
    for k = [1:s_y]
        s_out(k) = map(find(divided(:,k)== min_div(k)));
    end
end