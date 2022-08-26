function bf = brown_func(X)

% Function for computing the gradient vector in a specified point of Brown
% function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% bgrad: n by 1 vector representing the gradient in X

f = @(x, k) (x(k-1).^2).^(x(k).^2 + 1) + (x(k).^2).^(x(k-1).^2 + 1);

n = length(X);
bf = 0;

for i = 1:n/2
    j = 2*i;
    bf = bf + f(X, j);
end

end