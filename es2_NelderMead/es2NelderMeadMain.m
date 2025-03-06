% Exercise 2, Nelder Mead method
% In this script the Nelder Mead method is applied to minimize the Rosenbrock
% function in two dimensions with two different starting point. 

clc
clear all
close all

dim = 2;  
initial_point1 = [1.2,1.2];
initial_point2 = [-1,1.2];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % Rosenbrock function
x_opt = [1,1];

% Tuned parameters for point 1
rho1 = 0.7;
sigma1 = 0.1;
gamma1 = 0.2;
chi1 = 1.5;

% Tuned parameters for point 2
rho2 = 0.7;
sigma2 = 0.1;
gamma2 = 0.6;
chi2 = 2;

% tol for first stopping criterion
tol_simplex = 1e-07;  
% tol for second stopping criterion
tol_varf = 1e-15;  
% maximum number of iteration
kmax = 10000;

% Create initial simplex
[simplex_initial1, flag1] = NelderMead_simplex(dim, initial_point1);
[simplex_initial2, flag2] = NelderMead_simplex(dim, initial_point2);

% Nelder mead output:
% [# iteration, last simplex, flag about convergence, type of move through
% iteration (expansion, contraction, shrinking), max distance between
% barycenter and vertices through iterations]
tic;
[k1, simplex1, x_1, flag1, type_of_move1, size_vec_1]  = nelder_mead(f, simplex_initial1, kmax, rho1, chi1, gamma1, dim, sigma1, tol_simplex, tol_varf);
time1 = toc;
tic;
[k2, simplex2, x_2, flag2, type_of_move2, size_vec_2]  = nelder_mead(f, simplex_initial2, kmax, rho2, chi2, gamma2, dim, sigma2, tol_simplex, tol_varf);
time2 = toc;

% Output
disp("Regarding point of convergence...")
disp("convergence for [1.2,1.2] ?")
disp(flag1)
disp("convergence for [-1,1.2] ?")
disp(flag2)
disp("[1.2,1.2] initial point, convergence point")
disp(x_1(end,:))
disp("distance from the optimum:")
disp(norm(x_1(end,:) - x_opt))
disp("[-1,1.2] initial point, convergence point")
disp(x_2(end,:))
disp("distance from the optimum:")
disp(norm(x_2(end,:) - x_opt))
disp("  ")
disp("Regarding speed of convergence...")
disp("[1.2,1.2] initial point, number iteration before convergence")
disp(k1)
disp("[-1,1.2] initial point, number iteration before convergence")
disp(k2)
disp("[1.2,1.2] initial point, computational costs")
disp(time1)
disp("[-1,1.2] initial point, computational costs")
disp(time2)
disp("size of last simplex for point [1.2,1.2] ")
disp(size_vec_1(end))
disp("size of last simplex for point [-1,1.2]")
disp(size_vec_2(end))

% Theorical rates:
vec_rate1 = compute_errorRatio(x_1, k1, x_opt);   
%vec_rate2 = compute_errorRatio(x_2, k2, x_opt);

% Experimental rates:
disp("Last 5 exponential rate of initial point [1.2,1.2]:")
vec_rate3 = compute_exp_rate(x_1, k1);  
fprintf('%g ', vec_rate3(end-5:end));
fprintf('\n'); 
disp("Last 5 exponential rate of initial rates [-1,1.2]")
vec_rate4 = compute_exp_rate(x_2, k2);
fprintf('%g ', vec_rate4(end-5:end));
fprintf('\n'); 

% Stagnation:
vec_increments1 = stagnation_func(x_1);
vec_increments2 = stagnation_func(x_2);

% Picture
f = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
x_interval = linspace(-2, 2, 500);  
y_interval = linspace(-1, 3, 500); 
nelderMead_picture2D(f, x_interval, y_interval, x_1, initial_point1, type_of_move1)
nelderMead_picture2D(f, x_interval, y_interval, x_2, initial_point2, type_of_move2)

figure;
plot(1:length(size_vec_1), size_vec_1)
xlabel('iter');
ylabel('max distance vertice-x_bar');
title('Size of the simplex through iterations');

figure;
plot(1:length(size_vec_2), size_vec_2)
xlabel('iter');
ylabel('max distance vertice-x_bar');
title('Size of the simplex through iterations');

figure;
plot(1:length(x_1), vecnorm(x_1 - x_opt, 2, 2),'LineWidth', 1.5,  'Color', 'b')
xlabel('iter','FontSize', 14);
ylabel('distance x-x_opt', 'FontSize', 14);
title('Distance from optimum through iterations','FontSize', 14);

figure;
plot(1:length(x_2), vecnorm(x_2 - x_opt, 2, 2),'LineWidth', 1.5,  'Color', 'b')
xlabel('iter','FontSize', 14);
ylabel('distance x-x_opt', 'FontSize', 14);
title('Distance from optimum through iterations','FontSize', 14);

