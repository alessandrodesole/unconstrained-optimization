function p76hess = problem_76_hess(X)

% Function for computing the Hessian matrix in a specified point of the
% function reported in problem 76 in Test Problems for Unconstrained Optimization

% INPUTS
% X: n by 1 vector representing a point in the n-dimensional space

% OUTPUTS
% p76hess: a n by n sparse tridiagonal matrix representing the Hessian
% matrix at point X

n = length(X);
max_nz = 3*n;
row = zeros(max_nz, 1);
column = zeros(max_nz, 1);
values = zeros(max_nz, 1);
nz = 1;
tmp = -0.2 * X(1);
row(nz) = n;
column(nz) = 1;
values(nz) = tmp;
nz = nz + 1;
row(nz) = 1;
column(nz) = n;
values(nz) = tmp;
nz = nz + 1;
row(nz) = 1;
column(nz) = 1;
values(nz) = -0.2 * X(n) + 0.06 * X(1) .^ 2 + 1;

for i = 2:n
    nz = nz + 1;
    row(nz) = i;
    column(nz) = i;
    values(nz) = -0.2 * X(i-1) + 0.06 * X(i).^ 2 + 1;
    tmp = -0.2 * X(i);
    nz = nz + 1;
    row(nz) = i;
    column(nz) = i - 1;
    values(nz) = tmp;
    nz = nz + 1;
    row(nz) = i-1;
    column(nz) = i;
    values(nz) = tmp;
end

p76hess = sparse(row, column, values, n, n);

end