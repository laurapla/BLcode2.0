function sigma = sigmav(tau_v, T_vl, Sa, df)

sigma = 1;

if Sa<0
    sigma = 4.0;
end
if tau_v>T_vl && tau_v<=2*T_vl
    if Sa<0
        sigma = 2.0;
    else
        sigma = 3.0;
    end
elseif df>0
    sigma = 4.0;
elseif Sa<0 && df>0
    sigma = 1.0;
end

end