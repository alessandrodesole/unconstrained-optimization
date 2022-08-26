function p81f = problem_81_func(X)

% Function for computing the gradient vector in a specified point of
% problem 81 function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% p81f: n by 1 vector representing the gradient in X

f_1 = @(x, k) x(k+1).^2 - 1
f_k = @(x,k) x(k).^2 + log(x(k+1)) - 1;

n = length(X);
p81f = 0;

for i = 1:n-1
    if i == 1
       p81f = p81f + f_1(X, i);
    else
       p81f = p81f + f_k(X, i);
    end
end

p81f = 0.5 * p81f;

end