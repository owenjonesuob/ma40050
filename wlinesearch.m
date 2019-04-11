function alpha = wlinesearch(x, f, df, s, theta_sd, theta_c)
% Backtracking line search algorithm which guarantees Wolfe conditions

% Check parameters are valid
if 0 >= theta_sd || theta_sd >= theta_c || theta_c >= 1
    error("Please specify 0 < theta_sd < theta_c < 1")
end
    
% Initialise alpha and auxiliary variables
alpha = 1;
a1 = 0;
a2 = 0;

% Define indicator functions for Wolfe convergence conditions
armijo = @(z, a) f(z + a*s) <= f(z) + theta_sd*a*(df(z)'*s);
curvature = @(z, a) df(z + a*s)'*s >= theta_c*(df(z)'*s);

% Adjust alpha until the conditions are satisfied; stop after a while
for k = 1:100

    if ~armijo(x, alpha)
        % Reduce alpha
        a2 = alpha;
        alpha = 0.5*(a1 + a2);
        
    elseif ~curvature(x, alpha)
        % Increase alpha
        a1 = alpha;
        if a2 == 0
            alpha = 2*a1;
        else
            alpha = 0.5*(a1 + a2);
        end
        
    else
        % We satisfied both conditions!
        break
        
    end
    
end

end
