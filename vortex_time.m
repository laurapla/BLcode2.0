function [tauv] = vortex_time(Cnprime,s_time,C_N1)

% Function that computes the force of the vortex that sheds from the
% leading edge as a function of the non-dimensional time
% Cv = force/circulation of the leading edge vortex
% tauv = non-dimensional time parameter counting the time from when the
% vortex sheds from the leading edge to when it reaches the trailing edge
% and sheds into the wake
% Cnc = circulatory normal force coefficient (attached flow)
% Cnprime = equivalent critical normal force coefficient
% s_time = non-dimensional time vector
% fprimeprime = unsteady trailing edge separation point
% T_vl = time constant for vortex traverse over chord [s]
% C_N1 = critical normal force coefficient

N = length(s_time);

% Non-dimensional time-conter (counts the passage of the vortex over the
% airfoil chord)
dtau = s_time(2);
tauv = zeros(1,N);

resta = abs(Cnprime)-C_N1;

% tauv starts counting when abs(Cnprime) exceeds C_N1, and goes back to 0
% when it exceeds the time of passage of the vortex over the airfoil
for i = 2:N
    if resta(i)>=0 && resta(i-1)<0
        tauv(i) = (s_time(i)-s_time(i-1))*(abs(Cnprime(i))-C_N1)/(abs(Cnprime(i))-abs(Cnprime(i-1)));
    elseif resta(i)>=0 && resta(i-1)>=0
        tauv(i) = tauv(i-1)+dtau;
    else
        tauv(i) = 0;
    end
end