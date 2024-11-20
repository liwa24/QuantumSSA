function theta = kepler_solver(type, e, tol, mu, varargin) 
% Kepler's inverse problem solver function (time or mean anomaly -> true
% anomaly) for a vector of n initial data
% the problem is solved using the built-in function 'fsolve'
%
% PROTOTYPE:
%   theta = kepler_solver(type, e, tol, mu, varargin)
%
% INPUT:
%   type          [string] Type of problem: 
%                   'from_t'  start from time instant
%                   'from_M'  start from mean anomaly
%   e             [n,1]     eccetricity of the orbit                []
%   tol           [1,1]     tolerance of the solver                 []
%   mu            [1,1]     gravitational parameter                 [km^3/s^2]
%   varargin      set of case-dependent input
%       case 'from_t'
%           t   [n,1] time vector                           [s]
%           a   [n,1] semi-major axis                       [km]
%           t0  [n,1] initial time instant (NOT mandatory)  [s]
%           th0 [n,1] initial true anomaly (NOT mandatory)  [rad]
%       case 'from_M'
%           M   [n,1] mean anomaly                          [rad]
%           th0 [n,1] initial true anomaly (NOT mandatory)  [rad]      
% 
%
% OUTPUT:
%   theta       [n,1]      true anomalies vector, solution of kep prob               [rad]
%
% VERSION:
% 01 - Walther Litteri, 11/20
% 02 - Walther Litteri, 01/24
% ---------------------------------------------------------------------

if sum(e<1) == length(e)
    if strcmp(type,'from_t')
        t = varargin{1};
        a = varargin{2};

        n = sqrt(mu/a^3);
        
        if length(varargin) < 3 
            M = n.*t;
            E0 = zeros(length(M),1);
        else
            t0 = varargin{3};
            th0 = varargin{4};
            M = n.*(t-t0);
            E0 = 2*atan(sqrt((1-e)./(1+e)).*tan(th0./2)); % inital eccentric anomaly
        end
        
    elseif strcmp(type,'from_M')
        M = varargin{1};

        if length(varargin) < 2 
            E0 = zeros(length(M),1);
        else
            th0 = varargin{2};
            E0 = 2*atan(sqrt((1-e)./(1+e)).*tan(th0./2));
        end

    end

    k = floor(M./(2*pi)); % it defines the number of revolutions
    Mbar = M - k*2*pi; %M parameter reduced to 1 revolution

%% solution of kepler problem wrt 1 revolution

Fsolve = @(E, M, e, E0)  E - e*sin(E)- M - (E0 - e*sin(E0));
options = optimset('TolX', tol);
Ebar = zeros(length(Mbar), 1);
theta = Ebar;

for i = 1: length(Mbar)
    Fsolve_i = @(E) Fsolve(E, Mbar(i), e, E0);
    Eiguess = Mbar(i) + (e*sin(Mbar(i)))/(1-sin(Mbar(i)+e) + sin(Mbar(i))); %cfr pag 71 lab 01
    Ebar(i) = fzero(Fsolve_i, Eiguess, options);  %fzero solutor has been used
    theta(i) = 2*atan(sqrt((1+e)/(1-e))*tan(0.5*Ebar(i))) ;
    if theta(i) < 0
        theta(i) = theta(i) + 2*pi;
    end
 theta(i) = theta(i) + 2*k(i)*pi;  %conversion to the number of revolutions k

end


elseif sum(e>1) == length(e)
     if strcmp(type,'from_t')
        t = varargin{1};
        a = varargin{2};
        
        n = sqrt(mu/abs(a)^3);

        if length(varargin) < 3 
            Mbar = n.*t;
            E0 = zeros(length(Mbar,1));  
        else
            t0 = varargin{3};
            th0 = varargin{4};
            Mbar = n.*(t-t0);
            E0 = 2*atanh(sqrt((e-1)/(e+1))*tan(th0/2)); 
        end
        
    elseif strcmp(type,'from_M')
        Mbar = varargin{1};

        if length(varargin) < 2 
            E0 = zeros(length(Mbar),1);
        else
            th0 = varargin{2};
            E0 = 2*atan(sqrt((1-e)./(1+e)).*tan(th0./2));
        end

     end

    

Fsolve = @(E, M, e, E0)  -E + e*sinh(E)- M + (E0 - e*sinh(E0));
options = optimset('TolX', tol);
Ebar = zeros(length(Mbar), 1);
theta = Ebar;
for i = 1: length(Mbar)
    Fsolve_i = @(E) Fsolve(E, Mbar(i), e, E0);
    Eiguess = Mbar(i) ; %+ (e*sin(Mbar(i)))/(1-sin(Mbar(i)+e) + sin(Mbar(i))); %cfr pag 71 lab 01
    Ebar(i) = fzero(Fsolve_i, Eiguess, options);  %fzero solutor has been used
    theta(i) = 2*atan(sqrt((1+e)/(e-1))*tanh(0.5*Ebar(i))) ;
    if theta(i) < 0
        theta(i) = theta(i) + 2*pi;
    end


end


else

    if strcmp(type,'from_t')
        t = varargin{1};
        a = varargin{2};
        

        n = sqrt(mu/abs(a)^3);

        if length(varargin) < 3 
            Mbar = n.*t;
            
        else
            t0 = varargin{3};
            Mbar = n.*(t-t0);
            
        end
        
    elseif strcmp(type,'from_M')
        Mbar = varargin{1};

     end
    
    
    Ebar = zeros(length(Mbar), 1);
    theta = Ebar;
    for i = 1: length(Mbar)
        Ebar(i) = (3*Mbar(i) + sqrt(1+(3*Mbar(i))^2))^(1/3);
        theta(i) = 2*atan(Ebar(i) - 1/Ebar(i)) ;
        if theta(i) < 0
            theta(i) = theta(i) + 2*pi;
        end
    end

end

