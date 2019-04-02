
function all_x = bfgs_ex(f, df, x0, B0, tol)
% Performs generalised steepest descent, using the BFGS method to update
% an approximate Hessian matrix and returning a matrix of all iterates

% Only permit a certain maximum number of iterations
max_iters = 1000;

% Ensure dimensions of inputs are correct
[m, n] = size(x0);
if n ~=  1
    error("x0 must be a n*1 vector")
elseif ~isequal(size(B0), [m, m])
    error("B0 must be square and have same dimension as x0")
end

% Prepare array to store iterates
all_x = NaN(m, max_iters+1);
all_x(:, 1) = x0;

% Set initial guesses
x = x0;
H = inv(B0);


for k = 1:max_iters
    
    % Check convergence
    if norm(df(x)) <= tol
        break
    end
    
    s = - H * df(x);
    
    % Implement SLOW exact line search: as long as f is decreasing,
    % increase alpha by a small amount
    %alpha = 0;
    %inc = 1e-6;
    %while f(x + alpha*s) > f(x + (alpha + inc)*s)
    %    alpha = alpha + inc;
    %end
    
    % Succinct version of the above: using MATLAB's `fzero()` to find
    %   inf{alpha >= 0 : p'(alpha) = 0}
    % where p(alpha) = f(x + alpha*s)
    alpha = fzero(@(a) s'*df(x + a*s), 0);
    
    % Calculate new iterate and useful vectors
    d = alpha .* s;
    y = x;
    x = x + d;
    y = df(x) - df(y);
    
    % Store new iterate
    all_x(:, k+1) = x;
    
    % Update H using SMW/BFGS method
    rho = 1 / (y'*d);
    H = (eye(m) - rho*(d*y')) * H * (eye(m) - rho*(y*d')) + rho*(d*d');
    
end

% Return all iterates
all_x = all_x(:, 1:k);

end