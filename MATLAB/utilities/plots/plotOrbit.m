function [X, Y, Z] = plotOrbit(kepE1, mu,Thin, Thf, stepTh)
% Plot the arc length of the orbit described by a set of orbital
% elements for a specific arc length.
% INPUT:
% kepEl [1x5] orbital elements [km,rad]
% mu [1x1] gravitational parameter [km^3/s^2]
% Thin initial true anomaly [rad]
% Thf final true anomaly [rad]
% stepTh [1x1] arc length step [rad]
%
% OUTPUT:
% X [1xn] X position [km]
% Y [1xn] Y position [km]
% Z [1xn] Z position [km]

if Thf<Thin; Thf = Thf + 2*pi; end  % control on Thin vs Thf
th = Thin: stepTh:Thf;
X = zeros(length(th), 1); Y = zeros(length(th), 1); Z = zeros(length(th), 1);
for k = 1: length(th)
    rth = kep2car([kepE1, th(k)], mu);
    X(k) = rth(1); Y(k) = rth(2); Z(k) = rth(3);
end
    