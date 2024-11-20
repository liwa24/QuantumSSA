function X_dot = eq_motion(t, X, mu)
% Function 'EQ_MOTION' is the input of ODE solver.
%
% Maria Anna Laino - Oct 2021

    %----------------------------------------------------------------------
    % DEFINING THE INPUT:
    x = X(1);
    y = X(2);
    z = X(3);
    v_x = X(4);
    v_y = X(5);
    v_z = X(6);
    
    % Position of spacecraft wrt primaries
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x-(1-mu))^2 + y^2 + z^2);
    
    %----------------------------------------------------------------------
    % ODE:
    x_dot  = v_x;
    y_dot  = v_y;
    z_dot  = v_z;
    x_ddot = x + 2*v_y - (1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3;
    y_ddot = y - 2*v_x - y*((1-mu)/r1^3 + mu/r2^3);
    z_ddot = - z*((1-mu)/r1^3 + mu/r2^3);
    
    %----------------------------------------------------------------------
    % DEFINING THE OUTPUT:
    X_dot = zeros(length(X), 1);
    
    X_dot(1) = x_dot;
    X_dot(2) = y_dot;
    X_dot(3) = z_dot;
    X_dot(4) = x_ddot;
    X_dot(5) = y_ddot;
    X_dot(6) = z_ddot;

end

