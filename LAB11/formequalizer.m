function [equalizer] = formequalizer(apha_b,tou,len)
    % 不支援重根計算
    % 等化器長度要自己調整
    syms z_inv n z;
    u(n) = heaviside(n);
    
    poly = []; 
    F = 0;
    for k = [1:length(apha_b)]      
        F = F + apha_b(k)*z_inv^(tou(k));
    end
    for k = [0:max(tou)]      
        if sum(k == tou)
            poly = [poly apha_b(k+1)];
        else
            poly = [poly 0];
        end
    end
    
    R = roots(poly)';
    F_ = partfrac(F^-1);
    F^-1;
    c_of_F = children(F_);
    
    equalizer = zeros(1,len);
    half_len = (len)/2;
    formula = 0;
    for k = [1:length(c_of_F)]    
        c_of_F{k};
        aa = subs(c_of_F{k},z_inv,z^-1);
        pole = double(solve(aa^-1==0));
        sol = iztrans(aa);
        
        if  (pole==1)
            warning("stable equalizer is not exist");
        elseif abs(pole) < 1 %ROC for stable
            formula = formula + sol*u(n);
            equalizer = equalizer + [zeros(1,half_len) double(subs(sol,n,[0:half_len-1]))];
        elseif abs(pole) > 1 %ROC for stable
            formula = formula - sol * subs(u(n),n,-n-1);
            equalizer = equalizer - [double(subs(sol,n,[-half_len:1:-1])) zeros(1,half_len)];
        end
    end
    formula
end