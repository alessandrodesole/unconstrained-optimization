function ergHess = extended_rosenbrock_hess(X)

% Function for computing the Hessian matrix in a specified point of the
% extended Rosenbrock

% INPUTS
% X: n by 1 vector representing a point in the n-dimensional space

% OUTPUTS
% ergHess: a n by n sparse tridiagonal matrix representing the Hessian
% matrix at point X

n = length(X);
max_nz = 2*n - 1;
row = zeros(max_nz, 1);
column = zeros(max_nz, 1);
values = zeros(max_nz, 1);
nz = 0;

for i = 1:(n-1)
    if mod(i,2)==1
        for j = 1:3
            nz = nz+1;
            switch j
                case 1
                    row(nz) = i;
                    column(nz) = i;
                    values(nz) = 600*X(i).^2 - 200 * X(i+1) + 1;
                case 2
                    row(nz) = i;
                    column(nz) = i+1;
                    values(nz) = -200 * X(i);
                case 3
                    row(nz) = i+1;
                    column(nz) = i;
                    values(nz) = -200 * X(i);
            end
            
        end
        
    else
        nz = nz+1;
        row(nz) = i;
        column(nz) = i;
        values(nz) = 100;
    end
end
nz = nz+1;
row(nz) = n;
column(nz) = n;
values(nz) = 100;
ergHess = sparse(row, column, values, n, n);
end