function [sigma] = sigma2(tau_v, T_vl, Sa, df)

sigma = 1;

if Sa<0
    sigma = 4.0;
elseif tau_v>T_vl && tau_v<=2*T_vl
    sigma = 3.0;
elseif df>0
    sigma = 4.0;
elseif tau_v>0 && tau_v<=T_vl && Sa<0
    sigma = 2.0;
elseif Sa<0 && sigma == 1.0
    sigma = 4.0;
elseif Sa<0 && df>0
    sigma = 1.0;
end

end