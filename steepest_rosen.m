% Set initial guesses
x0 = [-1.2; 1];
B0 = eye(2);

% Apply basic steepest descent algorithm
theta_sd = 0.1;
theta_c = 0.9;
tol = 1e-5;

x = bfgs_w(@f_rosen, @g_rosen, x0, B0, theta_sd, theta_c, tol);

% Set exact solution
xex = [1; 1];

% Visualise the solution path 
visual(@f_rosen, x, x0, xex);

n = size(x, 2);
fprintf("Number of iterations: %d\n", n);

e = vecnorm(x - xex);
fprintf("Final error: %e\n", e(n));


