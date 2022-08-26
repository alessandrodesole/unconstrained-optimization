function p81grad = problem_81_grad(X)

% Function for computing the gradient vector at a specified point of the
% function reported in problem 81 in Test Problems for Unconstrained Optimization

% INPUTS
% X: n by 1 vector represening a point in the n-dimensional space

% OUTPUTS
% p81grad: n by 1 vector which represent the gradient at X

n = length(X);
p81grad = zeros(n, 1);

p81grad1 = @(x, k) 4*x(k).^3 + 2 * x(k) * log(x(k+1)) - 4 * x(k);
p81gradn = @(x, k) (x(k-1).^2 + log(x(k)-1)) / x(k);
p81gradk = @(x, k) (x(k-1).^2 + log(x(k)) - 1) / x(k) + 2*x(k).^3 + 2 * x(k) * log(x(k+1)) - 2*x(k);

for i=1:n
    switch i
        case 1
            p81grad(i) = p81grad1(X, i);
        case n
            p81grad(i) = p81gradn(X, i);
        otherwise
            p81grad(i) = p81gradk(X, i);        
    end
end

end