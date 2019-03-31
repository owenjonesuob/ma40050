function out = f_quad(x)

A = [2, -1; -1, 10];
b = [-2; 1];

out = 0.5*x'*(A*x) + b'*x;

end