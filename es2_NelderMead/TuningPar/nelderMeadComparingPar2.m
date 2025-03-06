% Exercise 2: tuning parameters in Nelder Mead method for the Rosenbrock
%
% This script compares the performance of the Nelder-Mead optimization algorithm
% on the Rosenbrock function using different configurations of four parameters: 
% rho, sigma, gamma, and chi. Each parameter is tested with "l" different values.
% For each configuration of parameters (l x l x l x l in total), the number of iterations
% before convergence and the distance from the optimal solution are calculated.
% A weighted sum of the number of iterations and the convergence error is used 
% to determine the best configuration which is the one that minimize this weighted sum. 
% The idea is that the choice for the parameter is done taking into consideration both
% these two quantities: best parameter shoud minimize both of them.
% The preference for the weights depends on whether there is more interest in faster 
% convergence with less accuracy or slower convergence with higher precision.
%
% Outputs: 
% - Best configuration of parameters (rho, sigma, gamma, chi) that minimizes 
%   the weighted sum of iteration count and convergence error.

clc
clear all
close all

dim = 2;  
%initial_point = [-1.2,1];  
initial_point = [1.2,1.2];
x_opt = [1,1];
f = @(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2 ; % Rosenbrock function

% About parameters...
num_par = 4;
l = 9;  % number of possible value of each parameter to be tested
kmax = 10000;
rho_vec = [0.1, 0.3, 0.5, 0.7, 1, 1.3, 1.5, 1.7, 1.9];  % reflection
sigma_vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];  % shrinking
gamma_vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];   % contraction
chi_vec = [1.1, 1.3, 1.5, 1.7, 2, 2.1, 2.3, 2.5, 2.7];   % expansion
tol_simplex = 1e-07; 
tol_varf = 1e-15;  

% for choosing the best configuration
weight_k = 1*10^(-5);
weight_opt = 1 - weight_k;

% initial simplex
[simplex_initial, flag2] = NelderMead_simplex(dim, initial_point);


configuration_k = zeros(l,l,l,l); % store the number of iteration for each configuration of par
configuration_err_conv = zeros(l,l,l,l); % store the distance from optimal point for each configuration of par
for i_rho = 1:l
    for i_sig = 1:l
        for i_gam = 1:l
            for i_chi = 1:l
                [k, simplex,x, flag]  = nelder_mead(f, simplex_initial, kmax, rho_vec(i_rho), chi_vec(i_chi), gamma_vec(i_gam), dim, sigma_vec(i_sig), tol_simplex, tol_varf);
                configuration_k(i_rho, i_sig, i_gam, i_chi) = k;
                configuration_err_conv(i_rho, i_sig, i_gam, i_chi) = norm(x(end,:) - x_opt);
            end
        end
    end
end

% calculate the sum of number of iteration and distance from optimal point
% weighed.
configuration_qnt = (weight_k*configuration_k) + (weight_opt*configuration_err_conv);

% find the configuration which minimize the quantity 
[min_value, lin_index] = min(configuration_qnt(:));
[pos1, pos2, pos3, pos4] = ind2sub(size(configuration_qnt), lin_index);

disp("Best parameters:")
disp("rho")
disp(rho_vec(pos1))
disp("sigma")
disp(sigma_vec(pos2))
disp("gamma")
disp(gamma_vec(pos3))
disp("chi")
disp(chi_vec(pos4))


% Analysis of sensibility to parameters

% extreme values
max_k = max(configuration_k(:));
min_k = min(configuration_k(:));
max_err_conv = max(configuration_err_conv(:));
min_err_conv = min(configuration_err_conv(:));

% Number of iteration distribution varying parameter
figure;
histogram(configuration_k(:), 50);
title('Number of iteration distribution');
xlabel('Number of iterations');
ylabel('Frequency');
grid on;
