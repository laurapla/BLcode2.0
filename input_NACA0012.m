function [eta, C_Nalpha, alpha0, Cd0, C_M0, alpha1, dalpha1, S1, S2, K0, K1, K2, T_P, T_f, C_N1, T_v, T_vl, D_f] = input_NACA0012(M)

eta = 0.97; % Viscous effects factor

% Mach numbers for which there is data
M_data = [0.3 0.4 0.5 0.6 0.7 0.75 0.8];

C_Nalpha_data = [0.113 0.113 0.113 0.113 0.113 0.113 0.113]; % Normal force curve slope [1/deg]
alpha0_data = [0.17 0.17 0.17 0.17 0.17 0.17 0.17]; % Angle of zero lift [deg]
Cd0_data = [0 0 0 0 0 0 0]; % Drag coefficient at alpha = 0
C_M0_data = [-0.0037 -0.0037 -0.0037 -0.0037 -0.0037 -0.0037 -0.0037]; % Zero-lift moment coefficient

% Airfoil coefficients for unsteady separated flow modeling
alpha1_data = [15.25 12.5 10.5 8.5 5.6 3.5 0.7]; % Angle at which f=0.7 [deg]
dalpha1_data = [2.1 2.0 1.45 1.0 0.8 0.2 0.1]; % Dynamic offset of alpha1 [deg]
S1_data = [3.0 3.25 3.5 4.0 4.5 3.5 0.7]; % Coefficient that defines the stall characteristics [deg]
S2_data = [2.3 1.6 1.2 0.7 0.5 0.8 0.18]; % Coefficient that defines the stall characteristics [deg]
K0_data = [0.0025 0.006 0.02 0.038 0.03 0.001 -0.01]; % Aerodynamic center offset from the 1/4-chord
K1_data = [-0.135 -0.135 -0.125 -0.12 -0.09 -0.13 0.02]; % Direct effect on the center of pressure due to the growth of the separated flow region
K2_data = [0.04 0.05 0.04 0.04 0.15 -0.02 -0.01]; % Describes the shape of the moment break at stall
T_P_data = [1.7 1.8 2.0 2.5 3.0 3.3 4.3]; % Time constant for leading edge pressure response [s]
T_f_data = [3.0 2.5 2.2 2.0 2.0 2.0 2.0]; % Time constant for separation point [s]

% Airfoil coefficients for dynamic stall modeling
C_N1_data = [1.45 1.2 1.05 0.92 0.68 0.5 0.18]; % Critical normal force coefficient
T_v_data = [6.0 6.0 6.0 6.0 6.0 6.0 4.0]; % Time constant for vortex lift [s]
T_vl_data = [7.0 9.0 9.0 9.0 9.0 9.0 9.0]; % Time constant for vortex traverse over chord [s]
D_f_data = [8.0 7.75 6.2 6.0 5.9 5.5 4.0]; % Loss in chord force due to dynamic stall


% Interpolation
C_Nalpha = rad2deg(interp1(M_data,C_Nalpha_data,M)); % Normal force curve slope [1/rad]
alpha0 = deg2rad(interp1(M_data,alpha0_data,M)); % Angle of zero lift [rad]
Cd0 = interp1(M_data,Cd0_data,M); % Drag coefficient at alpha = 0
C_M0 = interp1(M_data,C_M0_data,M); % Zero-lift moment coefficient
alpha1 = deg2rad(interp1(M_data,alpha1_data,M)); % Angle at which f=0.7 [rad]
dalpha1 = deg2rad(interp1(M_data,dalpha1_data,M)); % Dynamic offset of alpha1 [rad]
S1 = deg2rad(interp1(M_data,S1_data,M)); % Coefficient that defines the stall characteristics [rad]
S2 = deg2rad(interp1(M_data,S2_data,M)); % Coefficient that defines the stall characteristics [rad]
K0 = interp1(M_data,K0_data,M); % Aerodynamic center offset from the 1/4-chord
K1 = interp1(M_data,K1_data,M); % Direct effect on the center of pressure due to the growth of the separated flow region
K2 = interp1(M_data,K2_data,M); % Describes the shape of the moment break at stall
T_P = interp1(M_data,T_P_data,M); % Time constant for leading edge pressure response [s]
T_f = interp1(M_data,T_f_data,M); % Time constant for separation point [s]
C_N1 = interp1(M_data,C_N1_data,M); % Critical normal force coefficient
T_v = interp1(M_data,T_v_data,M); % Time constant for vortex lift [s]
T_vl = interp1(M_data,T_vl_data,M); % Time constant for vortex traverse over chord [s]
D_f = interp1(M_data,D_f_data,M); % Loss in chord force due to dynamic stall

end