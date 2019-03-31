function g = g_rosen(x)
% The gradient of the two-dimensional Rosenbrock function 

b = 10;

g(1) = 2*(x(1) - 1) + 4*b*x(1)*(x(1)^2 - x(2));
g(2) = 2*b*(x(2) - x(1)^2);

g = g';

end

