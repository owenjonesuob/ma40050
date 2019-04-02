% Set initial guesses
x0 = [-1; 1];
B0 = eye(2);

% Apply basic steepest descent algorithm
tol = 1e-5;

x = bfgs_ex(@f_quad, @g_quad, x0, B0, tol);

% Set exact solution
xex = [1; 0];

% Visualise the solution path 
visual(@f_quad, x, x0, xex);

n = size(x, 2);
fprintf("Number of iterations: %d\n", n-1);

e = vecnorm(x - xex);
fprintf("Final error: %e\n", e(n));

