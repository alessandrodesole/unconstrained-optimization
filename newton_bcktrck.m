function [xk, fk, gradfk_norm, k, xseq, btseq, GMRESit] = ...
    newton_bcktrck(x0, f, gradf, Hessf, kmax, ...
    tolgrad, c1, rho, btmax)
%
% [xk, fk, gradfk_norm, k, xseq] = ...
%     newton_bcktrck(x0, f, gradf, Hessf, kmax, ...
%     tolgrad, c1, rho, btmax)
%
% Function that performs the newton optimization method, 
% implementing the backtracking strategy.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha0;
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy.
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% xseq = n-by-k matrix where the columns are the xk computed during the 
% iterations
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step.
%

% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) ...
    fk + c1 * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
btseq = zeros(1, kmax);
GMRESit = zeros(1, kmax);
fk = zeros(1, kmax+1);

xk = x0;
fk(1) = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = zeros(1, kmax+1);
gradfk_norm(1) = norm(gradfk);

while k < kmax && gradfk_norm(k+1) >= tolgrad
    % Compute the descent direction as solution of
    % Hessf(xk) p = - graf(xk)
    %Hessian = Hessf(xk);
    %[pk, flag, relres, iter] = gmres(Hessian, -gradfk);

    pk = -Hessf(xk)\gradfk;
    
    % Check for GMRES convergence
    if flag ~= 0
        warning("GMRES does not converge. Flag = %d", flag);
    end

    % Reset the value of alpha
    alpha = 1;
    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    bt = 0;

    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied
    while bt < btmax & fnew > farmijo(fk, alpha, gradfk, pk)
        % Reduce the value of alpha
        alpha = rho * alpha;
        % Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        % Increase the counter by one
        bt = bt + 1;
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    gradfk = gradf(xk);

    % Increase the step by one
    k = k + 1;

    gradfk_norm(k+1) = norm(gradfk);
    fk(k+1) = fnew;
    
    % Store current xk in xseq
    xseq(:, k) = xk;
    % Store bt iterations in btseq
    btseq(k) = bt;
    %GMRESit(k) = iter(2);
end

% "Cut" xseq and btseq to the correct size
xseq = xseq(:, 1:k);
btseq = btseq(1:k);
fk = fk(1:k+1);
gradfk_norm = gradfk_norm(1:k+1);
%GMRESit = GMRESit(1:k);

end