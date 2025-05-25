% Exercise 3, analysis of Modified Newton method 
% This script implements and tests the Modified Newton method
% simulations are performed to analyze the convergence, 
% computational cost, stagnation and empirical rates of convergence of the method in
% various dimensions (10^3,10^4,10^5) and multiple initial points (one suggested from 
% PDF and 10 points from the ipercube).
%% Compute graphics in order to compare results in different dimension %%
clc
clear all
load("output.mat")

% Rosenbrock chained %
create_bar_plot(mean(time_ros_1000), mean(time_ros_10000), mean(time_ros_100000), 'time', 'ros')
create_bar_plot(mean(diff_ros_1000), mean(diff_ros_10000), mean(diff_ros_100000), 'diff', 'ros')
create_bar_plot(mean(iter_ros_1000), mean(iter_ros_10000), mean(iter_ros_100000), 'iter', 'ros')
create_bar_plot(mean(final_rate_exact_deriv_ros_1000), mean(final_rate_exact_deriv_ros_10000), mean(final_rate_exact_deriv_ros_100000), 'rate', 'ros')

% Wood chained % 
create_bar_plot(mean(time_wood_1000), mean(time_wood_10000), mean(time_wood_100000), 'time', 'wood')
create_bar_plot(mean(diff_wood_1000), mean(diff_wood_10000), mean(diff_wood_100000), 'diff', 'wood')
create_bar_plot(mean(iter_wood_1000), mean(iter_wood_10000), mean(iter_wood_100000), 'iter', 'wood')
create_bar_plot(mean(final_rate_exact_deriv_wood_1000), mean(final_rate_exact_deriv_wood_10000), mean(final_rate_exact_deriv_wood_100000), 'rate', 'wood')

% Powell chained %
create_bar_plot(mean(time_pow_1000), mean(time_pow_10000), mean(time_pow_100000), 'time', 'pow')
create_bar_plot(mean(diff_pow_1000), mean(diff_pow_10000), mean(diff_pow_100000), 'diff', 'pow')
create_bar_plot(mean(iter_pow_1000), mean(iter_pow_10000), mean(iter_pow_100000), 'iter', 'pow')
create_bar_plot(mean(final_rate_exact_deriv_pow_1000), mean(final_rate_exact_deriv_pow_10000), mean(final_rate_exact_deriv_pow_100000), 'rate', 'pow')

%% Finite differences %%
% (1) fixed h %
clc
clear all
load("output.mat")

% Rosenbrock chained %
create_bar_plot(mean(time1_finite_diff_ros_1000,2), mean(time1_finite_diff_ros_10000,2), mean(time1_finite_diff_ros_100000,2), 'time', 'ros')
create_bar_plot(mean(norm_err_conv_ros_1000), mean(norm_err_conv_ros_10000), mean(norm_err_conv_ros_100000), 'diff', 'ros')
create_bar_plot(mean(iters_ros_1000), mean(iters_ros_10000), mean(iters_ros_100000), 'iter', 'ros')
create_bar_plot(final_rate_ros_1000, final_rate_ros_10000, final_rate_ros_100000, 'rate', 'ros')

% Wood chained % 
create_bar_plot(mean(time1_finite_diff_wood_1000,2), mean(time1_finite_diff_wood_10000,2), mean(time1_finite_diff_wood_100000,2), 'time', 'wood')
create_bar_plot(mean(norm_err_conv_wood_1000), mean(norm_err_conv_wood_10000), mean(norm_err_conv_wood_100000), 'diff', 'wood')
create_bar_plot(mean(iters_wood_1000), mean(iters_wood_10000), mean(iters_wood_100000,2), 'iter', 'wood')
create_bar_plot(final_rate_wood_1000, final_rate_wood_10000, final_rate_wood_100000, 'rate', 'wood')

% Powell chained %
create_bar_plot(mean(time1_finite_diff_pow_1000,2), mean(time1_finite_diff_pow_10000,2), mean(time1_finite_diff_pow_100000,2), 'time', 'pow')
create_bar_plot(mean(norm_err_conv_pow_1000), mean(norm_err_conv_pow_10000), mean(norm_err_conv_pow_100000), 'diff', 'pow')
create_bar_plot(mean(iters_pow_1000), mean(iters_pow_10000), mean(iters_pow_100000), 'iter', 'pow')
create_bar_plot(final_rate_pow_1000, final_rate_pow_10000, final_rate_pow_100000, 'rate', 'pow')

%% Finite differences %%
% (2) variable h %
clc
clear all
load("output.mat")

% Rosenbrock chained %
create_bar_plot(mean(time2_finite_diff_ros_1000,2), mean(time2_finite_diff_ros_10000,2), mean(time2_finite_diff_ros_100000,2), 'time', 'ros')
create_bar_plot(mean(norm_err_conv2_ros_1000), mean(norm_err_conv2_ros_10000), mean(norm_err_conv2_ros_100000), 'diff', 'ros')
create_bar_plot(mean(iters2_ros_1000), mean(iters2_ros_10000), mean(iters2_ros_100000), 'iter', 'ros')
create_bar_plot(final_rate2_ros_1000, final_rate2_ros_10000, final_rate2_ros_100000, 'rate', 'ros')


% Wood chained % 
create_bar_plot(mean(time2_finite_diff_wood_1000,2), mean(time2_finite_diff_wood_10000,2), mean(time2_finite_diff_wood_100000,2), 'time', 'wood')
create_bar_plot(mean(norm_err_conv2_wood_1000,2), mean(norm_err_conv2_wood_10000,2), mean(norm_err_conv2_wood_100000,2), 'diff', 'wood')
create_bar_plot(mean(iters2_wood_1000,2), mean(iters2_wood_10000,2), mean(iters2_wood_100000,2), 'iter', 'wood')
create_bar_plot(final_rate2_wood_1000, final_rate2_wood_10000, final_rate2_wood_100000, 'rate', 'wood')

% Powell chained %
create_bar_plot(mean(time2_finite_diff_pow_1000'), mean(time2_finite_diff_pow_10000'), mean(time2_finite_diff_pow_100000'), 'time', 'pow')
create_bar_plot(mean(norm_err_conv2_pow_1000), mean(norm_err_conv2_pow_10000), mean(norm_err_conv2_pow_100000), 'diff', 'pow')
create_bar_plot(mean(iters2_pow_1000), mean(iters2_pow_10000), mean(iters2_pow_100000), 'iter', 'pow')
create_bar_plot(final_rate2_pow_1000, final_rate2_pow_10000, final_rate2_pow_100000, 'rate', 'pow')
