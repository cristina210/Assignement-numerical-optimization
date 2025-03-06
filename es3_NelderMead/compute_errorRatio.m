% Function: compute_rate
% Computes the convergence rate of an iterative method given the convergence
% point. It outputs a vector of error ratios and plots the convergence 
% rates over iterations.
%
% Inputs:
% - x_bar: Matrix of solution estimates per iteration.
% - n_iter: Number of total iterations.
% - x_opt: Convergence point.
%
% Output:
% - vec_rate: Vector of convergence rates.

function [vec_rate] = compute_errorRatio(x_bar, n_iter, x_opt)
vec_rate = zeros(1,n_iter);
for i=1:(length(x_bar)-1)
    e_succ = norm(x_bar(i+1,:)-x_opt);
    e_k = norm(x_bar(i,:)-x_opt);
    ratio = e_succ/e_k;
    vec_rate(i) = ratio;
end
figure;
vec_rate = vec_rate(1:(length(x_bar)-1));
plot(1:length(vec_rate), vec_rate, 'LineWidth', 2)
title('Convergence rate','FontSize', 16)
xlabel('iteration', 'FontSize', 14)
ylabel('Error ratios', 'FontSize', 14)