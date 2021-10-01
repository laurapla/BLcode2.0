function [Cnp, Cmp, Ccp, alpha_E, Cn_i] = BL_attached(t_span,alpha,q,V,M,c,C_Nalpha,alpha0,x_ac)

% Function that computes the normal force coefficient, the moment
% coefficient and the drag force coefficient of the attached flow
% Cnp = normal force coefficient (vector)
% Cmp = moment coefficient (vector)
% Cdp = drag force coefficient (vector)
% alpha_E = effective angle of attack (Cnc/C_Nalpha) [rad]
% t_span = times at which the outputs will be computed (vector) [s]
% alpha = angle of attack (vector) [s]
% q = pitching rate (vector) [rad/s]
% V = free-stream velocity [m/s]
% M = Mach number
% c = airfoil chord [m]
% C_Nalpha = normal force curve slope [1/rad]
% x_ac = aerodynamic center normalized by the chord

% Time vector (to solve the dynamics)
time = t_span;

% Matrices of the system
[A,B,C,D] = attached_matrices(V,M,c,C_Nalpha,x_ac);

% Inputs of the system
u = [alpha; q];

% Initial conditions
xi = zeros(8,1);

% System of equations
[~,x] = ode45(@(t,x) dx_attached(t,x,A,B,alpha,q,time),t_span,xi);
x = x.';

% Output
y = C*x+D*u;
Cnp = y(1,:);
Cmp = y(2,:);

% Circulatory and non-circulatory components
alpha_E = effective_angle(sqrt(1-M^2),V,c,x(1,:),x(2,:)); % Effective angle [rad]
Cn_i = Cnp-C_Nalpha*(alpha_E-alpha0);

% Drag force
Ccp = C_Nalpha*sin(alpha_E-alpha0).^2; % Chord force coefficient

end