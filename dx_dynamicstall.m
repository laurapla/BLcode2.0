function dx = dx_dynamicstall(t,x,T_v,dCv,time)

% Function the gives the differential equation that accounts for the
% effects of the dynamic stall
% T_v = time cosntant for vortex lift [s]
% dCv = rate of change of circulation
% time = time instants at which dCv was computed

DCv = interp1(time,dCv,t);

dx = -x/T_v+DCv/T_v;

end