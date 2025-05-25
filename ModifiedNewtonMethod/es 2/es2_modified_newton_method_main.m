clc
clear all
close all

% set the seed to obtain the same results while running the program
seed_value = min(343341, 343428);
rng(seed_value);

dim = 2; 
initial_point1 = [1.2; 1.2];
initial_point2 = [-1.2; 1];
initial_points = [initial_point1, initial_point2];
x_opt = ones(2,1);

% function definition
f = @(x) 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;

% gradient of the function
gradf = @(x,h) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2; 
                    200*(x(2) - x(1)^2)];

% Hessian matrix of the function
hessf = @(x,h) [ 1200*x(1)^2 - 400*x(2) + 2, -400*x(1);
                    -400*x(1), 200];

% Implementation of centered finite differences approximation for gradf and hessf
gradf_finite_diff_rosenbrock = @(x, h)[
    -400*x(1)*x(2) + 400*x(1)^3 + 400*h(1).^2*x(1) - 2 + 2*x(1);
    200*x(2) - 200*x(1)^2
    ];

% in this way I obtain the hessian matrix for the case of h costante and
% variable, in the 1^st case I pass to the function h 2x1 with equal
% components
hess_finite_diff_rosenbrock = @(x, h)[
    1200*x(1)^2 - 400*x(2) + 2 + 200*h(1)^2,  -400*x(1) - 200*h(1);
    -400*x(1) - 200*h(1)               , 200
    ];

% parameters used by the function modified_netwon_method
rho = 0.5;
tolgrad = 10^(-7);
c1 = 10^(-4); 
btmax = 40;
kmax = 500;

flag_h = 0; % h not used, compute the exact derivatives
h = 0;

time = zeros(2,1);
xseq_tot = cell(1,2);
iter = zeros(2,1);
failure_tot = zeros(2,1);

% compute the experimental rate of convergence
vec_rate = zeros(1,2);
vec_rate_pre = zeros(1,2);

for i=1:2
    tic;
    [xk, fk, gradfk_norm, k, ...
        xseq, btseq, failure] = modified_newton_method(initial_points(:,i), f, gradf, hessf, kmax, ...
                                                                tolgrad, c1, btmax, dim, rho, h, flag_h);
    time(i) = toc;
    iter(i) = k;
    xseq_tot{i} = xseq;
    vec_rate(i) = compute_exp_rate_conv_multi(xseq(:,end-4:end));
    failure_tot(i) = failure;
end


% Picture
f_picture = @(x, y) 100*(y - x.^2).^2 + (1 - x).^2;
x_interval = linspace(-2, 2, 500);  
y_interval = linspace(-1, 3, 500);  
%modNewtonMethod_picture2D(f_picture, x_interval, y_interval, xseq_tot{1}, xseq_tot{2}, initial_point1, initial_point2)


% testing different parameters for bactracking strategy with rho ≠ 0.5 and
% c1 ≠ 10^-4 and starting from initial_point1
rho_vec = [0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9];
c1_vec = [10^(-7), 10^(-5), 10^(-4), 10^(-3), 0.01, 0.1];

time = zeros(length(rho_vec), length(c1_vec));
test_diff_c1_rho = zeros(length(rho_vec), length(c1_vec));
failure_tot = zeros(length(rho_vec), length(c1_vec));

time1 = zeros(length(rho_vec), length(c1_vec));
test_diff_c1_rho1 = zeros(length(rho_vec), length(c1_vec));
failure_tot1 = zeros(length(rho_vec), length(c1_vec));
for i=1:length(rho_vec)
    for j=1:length(c1_vec)
        tic;
        [xk, fk, gradfk_norm, k, ...
        xseq, btseq, failure] = modified_newton_method(initial_point1, f, gradf, hessf, kmax,...
                                                                tolgrad, c1, btmax, dim, rho_vec(i), h, flag_h);
        time(i,j) = toc;
        test_diff_c1_rho(i, j) = k;
        failure_tot(i,j) = failure;

        tic;
        [xk1, fk1, gradfk_norm1, k1, ...
        xseq1, btseq1, failure1] = modified_newton_method(initial_point2, f, gradf, hessf, kmax,...
                                                                tolgrad, c1, btmax, dim, rho_vec(i), h, flag_h);
        time1(i,j) = toc;
        test_diff_c1_rho1(i, j) = k1;
        failure_tot1(i,j) = failure1;
    end
end

% Find which rho and c1 values lead to the fastest convergence
min_value = min(test_diff_c1_rho(:));
[index_rho, index_c1] = find(test_diff_c1_rho == min_value);
min_rho = rho_vec(index_rho(1));
min_c1 = c1_vec(index_c1(1));
create_picture(test_diff_c1_rho, rho_vec, c1_vec, 1, 'c1', 'iter')
create_picture(test_diff_c1_rho1, rho_vec, c1_vec, 2, 'c1', 'iter')

% Finite differences


% parameter used for finite differences
h_vec = 10.^[-2, -4, -6, -8, -10, -12];


% Testing different values of h solving pcg without preconditioning
% starting from both initial_point1 and initial_point2
flag_h = 0; % we do not use a specific increment for each value x_i

time_finite_diff = zeros(length(rho_vec),length(h_vec));
norm_err_conv = zeros(length(rho_vec), length(h_vec));
iter = zeros(length(rho_vec), length(h_vec));
failure_tot = zeros(length(rho_vec),length(h_vec));
vec_rate = zeros(length(rho_vec),length(h_vec));

time_finite_diff2 = zeros(length(rho_vec),length(h_vec));
norm_err_conv2 = zeros(length(rho_vec), length(h_vec));
iter2 = zeros(length(rho_vec), length(h_vec));
failure_tot2 = zeros(length(rho_vec),length(h_vec));
vec_rate2 = zeros(length(rho_vec),length(h_vec));

% Testing the finite differences applying a specific h_i when
% differentiating with respect to the variable x_i if flag_h == 1
flag_h2 = 1;

time_finite_diff_h_var = zeros(length(rho_vec), length(h_vec));
norm_err_conv_h_var = zeros(length(rho_vec), length(h_vec));
iter_h_var = zeros(length(rho_vec), length(h_vec));
failure_tot_h_var = zeros(length(rho_vec),length(h_vec));
vec_rate_h_var = zeros(length(rho_vec),length(h_vec));

time_finite_diff_h_var2 = zeros(length(rho_vec), length(h_vec));
norm_err_conv_h_var2 = zeros(length(rho_vec), length(h_vec));
iter_h_var2 = zeros(length(rho_vec), length(h_vec));
failure_tot_h_var2 = zeros(length(rho_vec),length(h_vec));
vec_rate_h_var2 = zeros(length(rho_vec),length(h_vec));

tolgrad = 10^(-7);
c1 = 10^(-4); 
btmax = 40;
for i=1:length(rho_vec)
    for j=1:length(h_vec)
        tic;
        [xk, fk, gradfk_norm, k, ...
        xseq, btseq, failure] = modified_newton_method(initial_point1, f, gradf_finite_diff_rosenbrock, ...
                                                            hess_finite_diff_rosenbrock, kmax,...
                                                            tolgrad, c1, btmax, dim, rho, h_vec(j), flag_h);
        time_finite_diff(i,j) = toc;
        iter(i,j) = k;
        failure_tot(i,j) = failure;
        norm_err_conv(i,j) = norm(x_opt - xk,2);
        vec_rate(i,j) = compute_exp_rate_conv_multi(xseq);

        tic;
        [xk2, fk2, gradfk_norm2, k2, ...
        xseq2, btseq2, failure2] = modified_newton_method(initial_point2, f, gradf_finite_diff_rosenbrock, ...
                                                            hess_finite_diff_rosenbrock, kmax,...
                                                            tolgrad, c1, btmax, dim, rho, h_vec(j), flag_h);
        time_finite_diff2(i,j) = toc;
        iter2(i,j) = k2;
        failure_tot2(i,j) = failure2;
        norm_err_conv2(i,j) = norm(x_opt - xk2,2);
        vec_rate2(i,j) = compute_exp_rate_conv_multi(xseq2);

        tic;
        [xk_h_var, fk_h_var, gradfk_norm_h_var, k_h_var, ...
        xseq_h_var, btseq_h_var, failure_h_var] = modified_newton_method(initial_point1, f, gradf_finite_diff_rosenbrock, ...
                                                            hess_finite_diff_rosenbrock, kmax,...
                                                            tolgrad, c1, btmax, dim, rho, h_vec(j), flag_h2);
        time_finite_diff_h_var(i,j) = toc;
        iter_h_var(i,j) = k_h_var;
        failure_tot_h_var(i,j) = failure_h_var;
        norm_err_conv_h_var(i,j) = norm(x_opt - xk_h_var,2);
        vec_rate_h_var(i,j) = compute_exp_rate_conv_multi(xseq_h_var);

        tic;
        [xk_h_var2, fk_h_var2, gradfk_norm_h_var2, k_h_var2, ...
        xseq_h_var2, btseq_h_var2, failure_h_var2] = modified_newton_method(initial_point2, f, gradf_finite_diff_rosenbrock, ...
                                                            hess_finite_diff_rosenbrock, kmax,...
                                                            tolgrad, c1, btmax, dim, rho, h_vec(j), flag_h2);
        time_finite_diff_h_var2(i,j) = toc;
        iter_h_var2(i,j) = k_h_var;
        failure_tot_h_var2(i,j) = failure_h_var;
        norm_err_conv_h_var2(i,j) = norm(x_opt - xk_h_var,2);
        vec_rate_h_var2(i,j) = compute_exp_rate_conv_multi(xseq_h_var2);
    end
end

create_picture(iter, rho_vec, h_vec, 1, 'h', 'iter')
create_picture(iter_h_var, rho_vec, h_vec, 1, 'h\_var', 'iter')

