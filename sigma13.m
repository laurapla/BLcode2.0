function [sigma1, sigma3] = sigma13(Cnprime, C_N1, Sa, df, fprimeprime, fr, tauv, T_vl)

sigma1 = 1;
sigma3 = 5;

if abs(Cnprime)<C_N1
    if df<=0
        sigma1 = 1;
        sigma3 = 1;
    else
        sigma1 = 0.5;
        sigma3 = 5;
    end
elseif abs(Cnprime)>=C_N1
    if df<=0
        sigma1 = 1.75;
        sigma3 = 1.75;
    else
        sigma1 = 1.0;
        sigma3 = 5;
        if tauv>0 && tauv<=T_vl
            sigma1 = 0.25;
            if Sa>0
                sigma1 = 0.75;
            end
        end
    end
end

if abs(Cnprime)>C_N1 && df<=0
    if Sa<0 || fprimeprime<=0.7 || fr<=0.7
        sigma1 = 2.0;
        sigma3 = 2.0;
    end
end
    
    

end