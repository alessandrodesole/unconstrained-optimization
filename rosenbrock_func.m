function rf = rosenbrock_func(X)

% Function for computing the gradient vector in a specified point of
% Rosenbrock
% function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% bgrad: n by 1 vector representing the gradient in X

f = @(x, k) 100*(x(k,:)-x(k-1,:).^2).^2+(1-x(k-1,:)).^2;

n = length(X);
rf = 0;

for i = 2:n
    rf = rf + f(X, i);
end
end