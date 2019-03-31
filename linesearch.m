function alpha = linesearch(x, f, df, s, theta)
% Backtracking line search algorithm

% Initialise alpha and set function increment for sufficient descent 
% condition
alpha = 1;
delta = theta*df(x)'*s;

% Half alpha until the sufficient descent condition is satisfied; stop when
% alpha = 2^{-30} < 10^{-9}
for k = 1:30

    if f(x + alpha*s) - f(x) < alpha*delta
        break
    end
    
    alpha = alpha/2;
    
end

end
