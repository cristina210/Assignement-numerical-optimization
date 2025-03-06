% Function: comparePar2
% This function compares the performance of the Nelder-Mead optimization algorithm
% on the Rosenbrock function using different configurations of four parameters: 
% rho, sigma, gamma, and chi. Each parameter is tested with "l" different values.
% For each configuration of parameters (l x l x l x l in total), the number of iterations
% before convergence and the distance from the optimal solution are calculated.
% A weighted sum of the number of iterations and the convergence error is used 
% to determine the best configuration which is the one that minimize this weighted sum. 
%
% Inputs:
% - dim: the dimensiona of the problem 
% - f: the objective function to minimize
% - initial_point: The initial point for the Nelder-Mead algorithm
% - x_opt: The known optimal solution 
% - l: The number of different values to test for each parameter (rho, sigma, gamma, chi)
% - rho_vec, sigma_vec, gamma_vec, chi_vec: vectors containing the possible values for rho, sigma, gamma, and chi, respectively
%
% Outputs:
% pos1, pos2, pos3, pos4 - Indices of the best parameter configuration that minimizes the
% weighted sum of the number of iterations and the distance from the optimal solution.
% - pos1 corresponds to the index of the best rho value.
% - pos2 corresponds to the index of the best sigma value.
% - pos3 corresponds to the index of the best gamma value.
% - pos4 corresponds to the index of the best chi value.
%
% This function is used in "Tuning_parameter"

function [pos1, pos2, pos3, pos4] = comparePar2(dim, f, initial_point,x_opt, l,rho_vec,sigma_vec, gamma_vec,chi_vec)

% for Nelder Mead
kmax = 5000;
tol_simplex = 1e-07; 
tol_varf = 1e-07; 

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

