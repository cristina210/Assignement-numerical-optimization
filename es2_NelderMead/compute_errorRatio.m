% Function: compute_rate
% Computes the error ratio of an iterative method given the convergence
% point. It outputs a vector of error ratios and plots it.
%
% Inputs:
% - x: Matrix of solution estimates per iteration.
% - n_iter: Number of total iterations.
% - x_opt: Convergence point.
%
% Output:
% - vec_rate: Vector of convergence rates.

function [vec_rate] = compute_errorRatio(x, n_iter, x_opt)
n = length(x);
vec_rate = zeros(1, n-1);
for i = 1:(n-1)
    e_succ = norm(x(i+1,:) - x_opt);
    e_k = norm(x(i,:) - x_opt);
    if e_k ~= 0
        vec_rate(i) = e_succ / e_k;
    end
end
figure;
plot(1:length(vec_rate), vec_rate, 'LineWidth', 2);
title('Error ratio','FontSize', 16)
xlabel('iteration', 'FontSize', 14)
ylabel('Error ratios', 'FontSize', 14)