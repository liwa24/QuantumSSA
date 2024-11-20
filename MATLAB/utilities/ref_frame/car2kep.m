function [kepv] = car2kep(r, v, mu)
% conversione da descrizione cartesiana a parametri kepleriani
% INPUT:
% r [3x1] Position vector [km]
% v [3x1] Velocity vector [km/s]
% mu [1x1] Gravitational parameter [km^3/s^2]
%
% OUTPUT:
% a [1x1] Semi-major axis [km]
% e [1x1] Eccentricity [-]
% i [1x1] Inclination [rad]
% OM [1x1] RAAN [rad] - ascensione retta
% om [1x1] Pericentre anomaly [rad]
% th [1x1] True anomaly [rad]

rm = norm(r); vm= norm(v);  %definizione dei moduli del vettore posizione e velocità
E = 0.5*(vm^2)- mu/rm; a= -mu/(2*E); %calcolo del semiasse maggiore
h = cross(r,v); hm = norm(h); %calcolo del momento angolare e del modulo
ev = 1/mu*(cross(v,h)-mu*r./rm); e = norm(ev); %calcolo dell'eccentricità;
i = acos(1/hm*dot(h,[0,0,1])); %calcolo dell'inclinazione
N = cross([0,0,1],h); n = N/norm(N); % versore dell'asse dei nodi

if  isnan(n)
    OM = 0;
    n = [1; 0; 0];
else
    if n(2)>=0 %calcolo dell'ascensione retta
        OM=acos(n(1));
    else
        OM= (2*pi - acos(n(1)));
    end
end   
    if ev == zeros(3,1)
        om = 0;
    elseif  ev(3) > 0 %calcolo dell'anomalia del pericentro
        om = (acos(1/e*dot(n, ev)));
    else
        om = (2*pi - acos(1/e*dot(n, ev)));
    end 



if dot(v,r) >=0 %calcolo dell'anomalia vera 
    th = (acos(1/(rm*e)*dot(r,ev)));
else
    th = (2*pi - acos(1/(rm*e)*dot(r,ev)));
end
kepv = [a,e,i,OM,om,th];