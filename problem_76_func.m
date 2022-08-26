function p76f = problem_76_func(X)

% Function for computing the gradient vector in a specified point of
% problem 76 function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% p76f: n by 1 vector representing the gradient in X

f_k = @(x, k) x(k) - (x(k+1).^2)/10;
f_n = @(x,k) x(k) - (x(1).^2)/10;

n = length(X);
p76f = 0;

for i = 1:n
    if i < n
       p76f = p76f + f_k(X, i);
    else
        p76f = p76f + f_n(X, i);
    end
end

p76f = 0.5 * p76f;
end