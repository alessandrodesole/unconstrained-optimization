function rHess = rosenbrock_hess(X)

% Function for computing the Hessian matrix in a specified point of
% Rosenbrock
% function

% INPUTS
% X: n by 1 vector representing a point in n-dimensional space

% OUTPUTS
% rHess: n by 1 vector containing the Hessian matrix in X

%f = @(x, k) [1200*x(k-1, :).^2 - 400*x(k, :) + 2, -400*x(k-1, :); -400*x(k-1, :), 200];

fh11 = @(x, k) 1200*x(k).^2 - 400*x(k+1) + 2;
fh12 = @(x, k) -400*x(k);
fh21 = @(x, k) -400*x(k);
fh22 = @(x, k) 200;

n = length(X);
rHess = zeros(n, n);
 
for k = 1:2:n
    for i = k:k+1
        for j = k:k+1
            if mod(i, 2) == 1
                if mod(j, 2) == 1
                    rHess(i, j) = fh11(X, k);
                else
                    rHess(i, j) = fh12(X, k);
                end
            else
                if mod(j, 2) == 1
                    rHess(i, j) = fh21(X, k);
                else
                    rHess(i, j) = fh22(X, k);
                end
            end
         end
    end
end

end