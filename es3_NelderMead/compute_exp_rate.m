% Function: compute_exp_rate
% Computes the experimental rate of convergence
%
% INPUT:
%   - x: solutions for each iteration 
%   - n_iter: Total number of iterations.
%
% OUTPUT:
%   - vec_rate (vector 1 x k): Vector containing the computed experimental 
%     rates of convergence
%
%   The function calculates the experimental rate of convergence using the formula:
%       rate = log(e_succ / e_k) / log(e_k / e_prev)
%   where:
%       - e_succ is the distance between the current and next iteration.
%       - e_k is the distance between the current and previous iteration.
%       - e_prev is the distance between the previous and the one before it.
%

function [vec_rate] = compute_exp_rate(x_bar, n_iter)
vec_rate = zeros(1,n_iter);
for i=3:(length(x_bar)-1)
    e_succ = norm(x_bar(i+1,:)-x_bar(i,:));
    e_k = norm(x_bar(i,:)-x_bar(i-1,:));
    e_prev = norm(x_bar(i-1,:)-x_bar(i-2,:));
    num = log(e_succ/e_k);
    den = log(e_k/e_prev);
    vec_rate(i) = num/den;
end
vec_rate = vec_rate(3:(length(x_bar)-1));


