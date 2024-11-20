function [mu, LU,TU,VU,LPs] = constants_3BP(type)
switch type 
    case "EM"
        mu = 1.215e-2;
        LU = 389703;   % km  
        VU = LU/382981;     % km/s
        TU = 382981/(24*3600);     % days
        LPs = lagrange_points(mu);

       

end 
end

function [X] = lagrange_points(mu)
% This function computes the position of the Lagrangian points in
% normalized units

    
    % Lagrange point L1

        poly = [1, -(3-mu), (3-2*mu), -mu, 2*mu, -mu];
        rt = roots(poly);
        for k = 1:5
            if isreal(rt(k)) 
                d = rt(k); 
                X(1,1) = 1 - mu - d;
                X(2:3,1) = [0, 0];
            end
        end 
        
    % Lagrange point L2    
   
        poly = [1, (3-mu), (3-2*mu), -mu, -2*mu, -mu ];
        rt = roots(poly);
        for k = 1:5
            if isreal(rt(k)) 
                d = rt(k); 
                X(1,2) = 1 - mu + d;
                X(2:3,2) = [0, 0];
            end
        end
        
    % Lagrange point L3    
    
        poly = [1, (2+mu), (1+2*mu), -(1-mu), -2*(1-mu), -(1-mu)];
        rt = roots(poly);
        for k = 1:5
            if isreal(rt(k)) 
                d = rt(k); 
                X(1,3) = - mu - d;
                X(2:3,3) = [0, 0];
            end
        end 
        
    % Lagrange point L4    
    
        X(1,4) = 0.5 - mu;
        X(2,4) = sqrt(3)/2;
        X(3,4) = 0;
        
        
    % Lagrange point L5    
    
        X(1,5) = 0.5 - mu;
        X(2,5) = -sqrt(3)/2;
        X(3,5) = 0;   


end