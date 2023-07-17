function [sigmafN, sigmafM] = sigmaf(Cnprime, C_N1, Sa, df, fprimeprime, fr, tauv, T_vl)

% Function that accounts for the modification of the constant T_f when
% calculating the unsteady separation points for the normal force
% coefficient and the aerodynamic moment
% sigmafN = modifies T_f when calculating the normal force coefficient
% sigmafM = modifies T_f when calculating the aerodynamic moment
% Cnprime = equivalent critical normal force coefficient
% C_N1 = critical normal force coefficient
% Sa = difference in the angle of attack between this time step and the
% previous time step [rad]
% df = difference in the unsteady separation point normalized by the chord
% between this time step and the previous time step
% fprimeprime = unsteady separation point normalized by the chord
% fr = unsteady separation point normalized by the chord used in the
% calculation of the aerodynamic moment
% tauv = non-dimensional time parameter counting the time from when the
% vortex sheds from the leading edge to when it reaches the trailing edge
% T_vl = time constant for vortex traverse over chord

sigmafN = 1;
sigmafM = 5;

if abs(Cnprime)<C_N1
    if df<=0
        sigmafN = 1;
        sigmafM = 1;
    else
        sigmafN = 0.5;
        sigmafM = 5;
    end
elseif abs(Cnprime)>=C_N1
    if df<=0
        sigmafN = 1.75;
        sigmafM = 1.75;
    else
        sigmafN = 1.0;
        sigmafM = 5;
        if tauv>0 && tauv<=T_vl
            sigmafN = 0.25;
            if Sa>0
                sigmafN = 0.75;
            end
        end
    end
end

if abs(Cnprime)>C_N1 && df<=0
    if Sa<0 || fprimeprime<=0.7 || fr<=0.7
        sigmafN = 2.0;
        sigmafM = 2.0;
    end
end

end