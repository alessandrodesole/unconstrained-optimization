function erg = extended_rosenbrock_grad(X)

% Function for computing the gradient vector in a specified point of the
% extended Rosenbrock

% INPUTS
% X: n by 1 vector representing a point in the n-dimensional space

% OUTPUTS
% erg: n by 1 vector which represents the gradient at X

grad_odd = @(x, k) 200*x(k).^3 - 200*x(k)*x(k+1)+x(k)-1;
grad_even = @(x, k) 100*x(k) - 100*x(k-1).^2;

n = length(X);
erg = zeros(n, 1);

for i = 1:n
    if mod(i, 2) == 1
        erg(i) = grad_odd(X, i);
    else
        erg(i) = grad_even(X, i);
    end
end

end