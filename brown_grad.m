function bgrad = brown_grad(X)

% Function for computing the gradient vector in a specified point of Brown
% function

% INPUTS
% X: n by 1 vector representing a point in n-dimensional space

% OUTPUTS
% bgrad: n by 1 vector representing the gradient in X

%gradf = @(x, k) [2*x(k,:)*(x(k+1,:).^2 + 1)*(x(k,:).^2).^(x(k+1,:).^2) + 2*x(k,:)*log(x(k+1,:).^2)*(x(k+1,:).^2).^(x(k,:).^2+1); 2*x(k+1,:)*(x(k,:).^2 + 1)*(x(k+1,:).^2).^(x(k,:).^2) + 2*x(k+1,:)*log(x(k,:).^2)*(x(k,:).^2).^(x(k+1,:).^2+1)];

gradf1 = @(x, k) 2*x(k,:)*(x(k+1,:).^2 + 1)*(x(k,:).^2).^(x(k+1,:).^2) + 2*x(k,:)*log(x(k+1,:).^2)*(x(k+1,:).^2).^(x(k,:).^2+1);
gradf2 = @(x, k) 2*x(k+1,:)*(x(k,:).^2 + 1)*(x(k+1,:).^2).^(x(k,:).^2) + 2*x(k+1,:)*log(x(k,:).^2)*(x(k,:).^2).^(x(k+1,:).^2+1);

n = length(X);
bgrad = zeros(n, 1);

for i = 1:n
    if i < n
        k = i;
    else
        k = n-1;
    end
    if mod(i, 2) == 1
        bgrad(i) = gradf1(X, k);
    else
        bgrad(i) = gradf2(X, k);
    end
end
end