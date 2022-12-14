function p76grad = problem_76_grad(X)

% Function for computing the gradient vector at a specified point of the
% function reported in problem 76 in Test Problems for Unconstrained Optimization

% INPUTS
% X: n by 1 vector represening a point in the n-dimensional space

% OUTPUTS
% p76grad: n by 1 vector which represent the gradient at X

n = length(X);
p76grad = zeros(n, 1);

p76grad(1) = X(1) - 0.1 * X(2).^2 - 0.2*X(1)*X(n)+ 0.02 * X(1).^3;
p76grad(n) = -0.2*X(n-1) * X(n) + 0.02 * X(n).^3 + X(n) - 0.1 * X(1).^2;

for i=2:(n-1)
    p76grad(i) = -0.2 * X(i-1) * X(i) + 0.02 * X(i).^3 + X(i) - 0.1 * X(i+1).^2;
end

end
