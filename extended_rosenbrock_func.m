function erf = extended_rosenbrock_func(X)

% Function for computing the gradient vector in a specified point of
% extended Rosenbrock function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% erf: n by 1 vector representing the gradient in X

f_odd = @(x, k) (10*(x(k+1,:).^2)-x(k+2,:)).^2;
f_even = @(x,k) (x(k,:)-1).^2;

%f = @(x, k) 100*(x(k,:)-x(k-1,:).^2).^2+(1-x(k-1,:)).^2;

n = length(X);
erf = 0;

for i = 1:n-2
    if mod(i, 2) == 1
        erf = erf + f_odd(X, i);
    else
        erf = erf + f_even(X, i);
    end
end

erf = 0.5 * erf;
end