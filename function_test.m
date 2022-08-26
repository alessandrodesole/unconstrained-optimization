%% LOADING THE VARIABLES FOR THE TEST

clear
close all
clc

alpha0 = 1;
kmax = 10000;
c1 = 1e-4;
tolgrad = 1e-12;
rho = 0.5;
btmax = 50;
n = 2; % Space dimension (choose among 2, 1e3, 1e4)

% Inexact Newton method parameters

FDgrad = ''; % Exact Gradient
FDHess = ''; % Exact Hessian
fterms = ['l', 's', 'q']; % Forcing terms (choose among 'l', 's', 'q')
h = 1e-8;
pcg_maxit = 50;

% ---------------------------------------------------------------------- %

elapsed_time_in = 0;
elapsed_time_n = 0;
elapsed_time_sd = 0;

grad_norm_last_step_in = 0;
grad_norm_last_step_n = 0;
grad_norm_last_step_sd = 0;

fk_last_step_in = 0;
fk_last_step_n = 0;
fk_last_step_sd = 0;

pcg_it = 0;

format long

disp(['Dimension n: ' num2str(n, '%.0e')]);

%% Select the desired function to optimize

func_to_opt = 'Brown';

% Starting point (x0)

x0 = zeros(n, 1);

switch func_to_opt
    case 'Brown'
        for i = 1:n
            if mod(i, 2) == 1
                x0(i) = -1.0;
            else
                x0(i) = 1.0;
            end
        end
    otherwise
        for i = 1:n
            if mod(i, 2) == 1
                x0(i) = -1.2;
            else
                x0(i) = 1.0;
            end
        end
end


switch func_to_opt
    case 'Rosenbrock'
        if n == 2
            x0 = [-1.2; 1.0]; % or x0 = [1.2; 1.2];
        end
        f = @(x) rosenbrock_func(x);
        gradf = @(x) rosenbrock_grad(x);
        Hessf = @(x) rosenbrock_hess(x);
    case 'Brown'
        if n == 2
            x0 = [-1.0; 1.0];
        end
        f = @(x) brown_func(x);
        gradf = @(x) brown_grad(x);
        Hessf = @(x) brown_hess(x);
    case 'Extended_Rosenbrock'
        f = @(x) extended_rosenbrock_func(x);
        gradf = @(x) extended_rosenbrock_grad(x);
        Hessf = @(x) extended_rosenbrock_hess(x);
    case 'Problem_76'
        f = @(x) problem_76_func(x);
        gradf = @(x) problem_76_grad(x);
        Hessf = @(x) problem_76_hess(x);

     case 'Problem_81'
        f = @(x) problem_81_func(x);
        gradf = @(x) problem_81_grad(x);
end

% Initialize and name figures

r = 0;
figure(1) % figure showing function value trend w.r.t. no. iteration k
figure(2) % figure showing norm gradient trend w.r.t. num iteration k
colors = ['r', 'g', 'y']; % colors for plots

xk_in = zeros(n, 3);
fk_in = zeros(3, 1);
elapsed_time_in = zeros(3, 1);

for j = 1:length(fterms)
    
    ft = fterms(j);
    
    switch ft
        % linear convergence
        case 'l'
            ft_string = 'Inexact Newton linear';
            forcing_term = @(gradfk, k) 0.5;
        case 's'
            ft_string = 'Inexact Newton superlinear';
            forcing_term = @(gradfk, k) min(0.5, sqrt(norm(gradfk)));
        case 'q'
            ft_string = 'Inexact Newton quadratic';
            forcing_term = @(gradfk, k) min(0.5, norm(gradfk));
        otherwise % use linear convergence by default
            ft_string = 'Inexact Newton linear';
            forcing_term = @(gradfk, k) 0.5;
    end
      
    fprintf('---------- %s ------------', ft_string);
    % Start time
    tic;

    [xk, fk_in, gradfk_norm_in, k, xseq_in, btseq_in] = ...
        innewton_general(x0, f, gradf, Hessf, kmax, ...
        tolgrad, c1, rho, btmax, FDgrad, FDHess, h, forcing_term, pcg_maxit);
    
    % End time
    elapsed_time_in(j) = toc;
    
    xk_in(:, j)=xk;
    fk_in_all(j)=fk_in;
    gradfk_norm_in_all(j) = gradfk_norm_in;
   

    % Compute stats
    
    method_in = ft_string;
    grad_norm_last_step_in = gradfk_norm_in(end);
    fk_last_step_in = fk_in(end);
    r = j+1;

    figure(1)
    yyaxis left
    semilogy(fk_in_all, 'LineWidth', 2, 'Color', colors(j), 'DisplayName', ft_string)
    hold on
    figure(2)
    yyaxis left
    semilogy(gradfk_norm_in, 'LineWidth', 2, 'Color', colors(j), 'DisplayName', ft_string)
    hold on

end

for t = 1:3
    
    disp('**** INEXACT NEWTON : START *****')

    disp(['xk: ', mat2str(xk_in(:, t)), ' (actual minimum: [0; 0]);'])
    disp(['f(xk): ', mat2str(fk_in(end)), ' (actual min. value: 0);'])
    disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
    disp(['Elapsed time: ', num2str(elapsed_time_in(t)), ';'])
    disp('************************************')
    
    disp('**** INEXACT NEWTON : FINISHED *****')
    disp('**** INEXACT NEWTON : RESULTS *****')
end

% Newton method with backtracking

disp('**** NEWTON: START *****')

tic;

[xk, fk_n, gradfk_norm_n, k, xseq_n, btseq_n, GMRESit] = ...
newton_bcktrck(x0, f, gradf, Hessf, kmax, ...
tolgrad, c1, rho, btmax);

elapsed_time_n = toc;

disp('**** NEWTON: FINISHED *****')
disp('**** NEWTON: RESULTS *****')

disp(['xk: ', mat2str(xk), ' (actual minimum: [0; 0]);'])
disp(['f(xk): ', mat2str(fk_n(end)), ' (actual min. value: 0);'])
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
disp(['Elapsed time: ', num2str(elapsed_time_n), ';'])
disp('************************************')

% Compute stats

method_n = 'Newton method';
grad_norm_last_step_n = gradfk_norm_n(end);
fk_last_step_n = fk_n(end);
r = r+1;

% Plotting graph

figure(1)
yyaxis left
semilogy(fk_n, 'LineWidth', 2, 'Color', 'b', 'DisplayName', method_n)
hold off
yyaxis left
ylabel('Value of the function')
xlabel('K iteration')
grid on
%xticks(1:M)
legend('Location', 'southwest')
title(['Function value Trend comparing all methods for ', num2str(n), ' dimensions'], 'FontSize', 10);


% Plot of Norm Gradient

figure(2)
yyaxis left
semilogy(gradfk_norm_n, 'LineWidth', 2, 'Color', 'b', 'DisplayName', method_n)
grid on
hold off
ylabel('Norm-Gradient')
xlabel('K iteration')
grid on
%xticks(1:M)
legend('Location', 'southwest')
figure(2)
title(['Norm Gradient Trend comparing all methods for ', num2str(n), ' dimensions'], 'FontSize', 10);

% ************************* Steepest Descent Method with backtracking
% ******************** %

disp('**** STEEPEST DESCENT: START *****')

% Start time
tic;

[xk, fk_sd, gradfk_norm_sd, k, xseq_sd, btseq_sd] = ...
steepest_desc_bcktrck(x0, f, gradf, alpha0, ...
kmax, tolgrad, c1, rho, btmax);

% End time
elapsed_time_sd = toc;

disp('**** STEEPEST DESCENT: FINISHED *****')
disp('**** STEEPEST DESCENT: RESULTS *****')
disp('************************************')


disp(['xk: ', mat2str(xk), ' (actual minimum: [0; 0]);'])
disp(['f(xk): ', mat2str(fk_sd(end)), ' (actual min. value: 0);'])
disp(['N. of Iterations: ', num2str(k),'/',num2str(kmax), ';'])
disp(['Elapsed time: ', num2str(elapsed_time_sd), ';'])
disp('************************************')


% Compute stats

method_sd = 'Steepest Descent method';
grad_norm_last_step_sd = gradfk_norm_sd(end);
fk_last_step_sd = fk_sd(end);

% Plot graphs

figure(3)
yyaxis left
semilogy(fk_sd, 'LineWidth', 2); % function value plot

hold on

semilogy(gradfk_norm_sd, '--', 'LineWidth', 1, 'Color', 'r'); %Norm Grad plot
ylabel('Norm of the gradient');
legend('Fk trend', 'GradFk trend', 'Location', 'northeast');
grid on
hold off
title(['Function value Trend Steepest Descent ', num2str(n), ' Dimensions'], 'Fontsize', 10);

% Print results

disp(['Results with N = ', num2str(n), ' Dimensions']);
disp(['%------------------', 'Results with N = ', num2str(n), ' Dimensions', ' ----------------%'])
disp([' Mehod  ', 'k    ', 'elapsed_time    ', 'grad_norm_last_step ', 'fk_last_step    '])


%% 2D PLOTS

if n == 2

    switch func_to_opt
        case 'Rosenbrock'

            % Creation of the meshgrid for the contour-plot
            [X, Y] = meshgrid(linspace(-6, 6, 500), linspace(-6, 25, 500));
            % Computation of the values of f for each point of the mesh
            Z = 100*(Y-X.^2).^2+(1-X).^2;
        otherwise
            % Creation of the meshgrid for the contour-plot
            [X, Y] = meshgrid(linspace(-4, 4, 500), linspace(-4, 4, 500));
            % Computation of the values of f for each point of the mesh
            Z = (X.^2).^(Y.^2 + 1) + (Y.^2).^(X.^2 + 1); % function 2
    end
    
    % Plots
    
    % Simple Plot
    fig1 = figure();
    % Contour plot with curve levels for each point in xseq
    [C1, ~] = contour(X, Y, Z);
    hold on
    % plot of the points in xseq
    plot([x0(1) xseq_sd(1, :)], [x0(2) xseq_sd(2, :)], '--*')
    hold off
    
    % More interesting Plot
    fig2 = figure();
    % Contour plot with curve levels for each point in xseq
    % ATTENTION: actually, the mesh [X, Y] is to coarse for plotting the last
    % level curves corresponding to the last point in xseq (check it zooming
    % the image).
    

    [C2, ~] = contour(X, Y, Z, fk_sd);
    hold on
    % plot of the points in xseq
    plot([x0(1) xseq_sd(1, :)], [x0(2) xseq_sd(2, :)], '--*')
    hold off
    
    % Barplot of btseq
    %fig3 = figure();
    %bar(btseq)
    
    % Much more interesting plot
    fig4 = figure();
    surf(X, Y, Z,'EdgeColor','none')
    hold on
    plot3([x0(1) xseq_sd(1, :)], [x0(2) xseq_sd(2, :)], [f(x0), f(xseq_sd)], 'r--*')
    hold off


end

