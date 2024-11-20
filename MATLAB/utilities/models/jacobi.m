function [E,J] = jacobi(X,mu)
% DESCRIPTION
% This funcition computes the Jacobi constant associated to the CR3BP
% PROTOTYPE
%
% INPUT
% X     [6,1]   State vector in the barycentric RF      [x;y;z;vx;vy;vz]
% mu    [1,1]   Gravitational parameter m2/(m1 + m2)    [#]
% OUTPUT
% E     [1,1]   Energy of the orbit                     [LU^2/TU^2]
% J     [1,1]   Jacobi constant J = -2*E                [LU^2/TU^2]
%
% DEPENDENCIES
%
% NOTES
%
% AUTHOR AND VERSION
%	Ver. 1 - W. Litteri - 03/2024

x = X(1);
y = X(2);
z = X(3);
xp = X(4);
yp = X(5);
zp = X(6);

mu1 = 1-mu;
mu2 = mu;

r1 = sqrt((x+mu2)^2 + y^2 + z^2);
r2 = sqrt((x-mu1)^2 + y^2 + z^2);

K = 0.5*(xp^2 + yp^2 + zp^2); %kinetic energy 
Ubar = -0.5*(x^2 + y^2) -mu1/r1 - mu2/r2 -0.5*mu1*mu2;

E = K + Ubar;
J = -2*E;

