function [A,B,C,D] = attached_matrices(V,M,c,C_Nalpha,x_ac)

% Function that gives the matrices of the system of equations of the
% attached flow
% V = free-stream velocity
% M = Mach number
% c = airfoil chord
% C_Nalpha = normal force curve slope
% x_ac = position of the aerodynamic center (normalized by the chord)

beta = sqrt(1-M^2);

A1 = 0.636; A2 = 1-A1; A3 = 1.5; A4 = -0.5; A5 = 1.0;
b1 = 0.339; b2 = 0.249; b3 = 0.25; b4 = 0.1; b5 = 0.5;
kNa = 0.75; kNq = 0.75; kMa = 0.8; kMq = 0.8;
% A1 = 0.3; A2 = 1-A1; A3 = 1.5; A4 = -0.5; A5 = 1.0;
% b1 = 0.14; b2 = 0.53; b3 = 0.25; b4 = 0.1; b5 = 0.5;
% kNa = 1; kNq = 1; kMa = 1; kMq = 1;

a = V/M; % Sound speed
T_I = c/a; % Basic non-circulatory time constant

% Constants
K_alpha = kNa/((1-M)+C_Nalpha*beta*M^2*(A1*b1+A2*b2));
K_q = kNq/(0.5*(1-M)+2*pi*beta*M^2*(A1*b1+A2*b2));
K_alphaM = kMa*(A3*b4+A4*b3)/(b3*b4*(1-M));
K_qM = 7*kMq/(15*(1-M)+3*pi*beta*M^2*A5*b5);

% Elements of matrix A
a11 = -b1*beta^2*2*V/c;
a22 = -b2*beta^2*2*V/c;
a33 = -1/(K_alpha*T_I);
a44 = -1/(K_q*T_I); % T_I or T_1?
a55 = -1/(b3*K_alphaM*T_I);
a66 = -1/(b4*K_alphaM*T_I);
a77 = -b5*beta^2*(2*V/c);
a88 = -1/(K_qM*T_I);

% Elements of matrix C
c11 = C_Nalpha*(2*V/c)*beta^2*A1*b1;
c12 = C_Nalpha*(2*V/c)*beta^2*A2*b2;
c13 = 4/M*a33;
c14 = 1/M*a44;
c21 = c11*(0.25-x_ac);
c22 = c12*(0.25-x_ac);
c25 = -1/M*A3*a55;
c26 = -1/M*A4*a66;
c27 = -C_Nalpha/16*A5*b5*beta^2*(2*V/c);
c28 = -7/(12*M)*a88;

% Matrices
A = diag([a11 a22 a33 a44 a55 a66 a77 a88]);
B = [1 1 1 0 1 1 0 0;
    0.5 0.5 0 1 0 0 1 1].';
C = [c11 c12 c13 c14 0 0 0 0;
    c21 c22 0 0 c25 c26 c27 c28];
D = [4/M 1/M;
    -(A3+A4)/M -7/(12*M)];

end