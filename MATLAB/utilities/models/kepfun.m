function dz = kepfun(~, z, mu)
%definition of function of the 2bodies equation to put inside the ode solver
% INPUTS
% t - time(s)
% z - vector of position and velocity components z = [x, y, z, xdot,ydot,zdot]
% mu - gravitational parameter [km^3/s^2]
%
% OUTPUTS
% dz - vector of z - derivative [xdot, ydot, zdot, xdotdot, ydotdot, zdotdot]
% this is required in the ode solutor
% -------------------------------------------------------------------

r = sqrt(z(1)^2 + z(2)^2 + z(3)^2); %distance modulus
dz = zeros(6,1);  % construction of dz by components
dz(1) = z(4);
dz(2) = z(5);
dz(3) = z(6);
dz(4) = -mu/r^3.*z(1);
dz(5) = -mu/r^3.*z(2);
dz(6) = -mu/r^3.*z(3);
end