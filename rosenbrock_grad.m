function rgrad = rosenbrock_grad(X)

% Function for computing the gradient vector in a specified point of
% Rosenbrock function

% INPUTS
% X: n by 1 vector representing a point in n-dimensuional space

% OUTPUTS
% rgrad: n by 1 vector representing the gradient in X

gradf1 = @(x, k) 2*x(k,:) - 400*x(k,:)*(-x(k,:).^2 + x(k+1,:))-2;
gradf2 = @(x, k) -200*x(k,:).^2 + 200*x(k+1,:);

n = length(X);
rgrad = zeros(n, 1);

for i = 1:n
    if i < n
        k = i;
    else
        k = n-1;
    end
    if mod(i, 2) == 1
        rgrad(i) = gradf1(X, k);
    else
        rgrad(i) = gradf2(X, k);
    end
end

end