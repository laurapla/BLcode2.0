function [C_Nalpha, alpha0, S1, S2] = airfoil_input(filename)


% Loading data from the file
A = readmatrix(filename);
alpha = deg2rad(A(:,1));
Cl = A(:,2);
Cd = A(:,3);

N = length(alpha);

% Pre-calculations
Cn = Cl.*cos(alpha)-Cd.*sin(alpha);
Cc = Cl.*sin(alpha)-Cd.*cos(alpha);


% Find static stall angle
i1 = find(abs(alpha)<0.01);
while Cl(i1)>=Cl(i1-1)
    alpha_ss = alpha(i1);
    i1 = i1+1;
end


% Calculate normal force slope
tolerance = 0.60;
i2 = find(abs(alpha+alpha_ss*tolerance)<0.01);
i3 = find(abs(alpha-alpha_ss*tolerance)<0.01);
y = Cn(i2:i3);
x = alpha(i2:i3);
X = [ones(length(x),1) x];
b = X\y;
C_Nalpha = b(2);
Cn0 = b(1);
alpha0 = -Cn0/C_Nalpha;


% Calculate the separation breaking point
f = zeros(N,1);
for i = 1:N
    f(i) = (2*sqrt(abs(Cn(i))/(abs(alpha(i)-alpha0)*C_Nalpha))-1)^2;
    if f(i)>=1
        f(i) = 1-1e-12;
    elseif f(i)<=0
        f(i) = 1e-12;
    end
end

i4 = find(abs(alpha)<0.01);
% ii1 = find(abs(f(i4:end)-0.7)<0.04,1);
% alpha1 = alpha(ii1+i4-1);
alpha1 = alpha_ss*0.87;



% Calculate S1 and S2 coefficients using a least squares regression
A = 1.0; B = 0.3; C = 0.04; D = 0.66;
i5 = find(abs(alpha-alpha1)<0.01);
alp1 = alpha(i4:i5);
f1 = f(i4:i5);
x1 = abs(alp1)-alpha1;
y1 = abs((f1-A)/(-B));
b1 = lsqr(x1,log(y1));
S1 = 1/b1;

i6 = find(abs(alpha-deg2rad(26))<0.02);
alp2 = alpha(i5:i6);
f2 = f(i5:i6);
x2 = alpha1-abs(alp2);
y2 = abs((f2-C)/D);
b2 = lsqr(x2,log(y2));
S2 = 1/b2;


end