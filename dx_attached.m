function dx = dx_attached(t,x,A,B,alpha,q,time)

% Function that gives the system of equations of the attached flow
% A = coefficient matrix that is multiplied to the state vector x
% B = coefficient matrix that is multiplied to the input vector u
% alpha = vector of angle of attack as a function of time [rad]
% q = vector of pitching rate as a function of time [rad/s]
% time = time steps at which alpha and q were computed [s]

u = zeros(2,1);
u(1) = interp1(time,alpha,t);
u(2) = interp1(time,q,t);

dx = A*x+B*u;
% dx(1,:) = A(1,1)*x(1)+u(1)+u(2)/2;
% dx(2,:) = A(2,2)*x(2)+u(1)+u(2)/2;
% dx(3,:) = A(3,3)*x(3)+u(1);
% dx(4,:) = A(4,4)*x(4)+u(2);
% dx(5,:) = A(5,5)*x(5)+u(1);
% dx(6,:) = A(6,6)*x(6)+u(1);
% dx(7,:) = A(7,7)*x(7)+u(2);
% dx(8,:) = A(8,8)*x(8)+u(2);

end