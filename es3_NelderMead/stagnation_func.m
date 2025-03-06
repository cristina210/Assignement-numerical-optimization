% Function: stagnation_func
% This function calculates and visualizes the norm of the differences
% between consecutive simplex's barycenters. It is intended to 
% study the presence of a "stagnation point" by plotting the changes in 
% successive iterations.
%
% Inputs:
% - x_bar: A matrix where each row represents an iteration point, and 
%   columns represent value of a barycenter's coordinates.
%
% Outputs:
% - vec_e_k: A row vector containing norms of the differences 
%   between consecutive simplex's barycenters.

function [vec_e_k] = stagnation_func(x_bar)
vec_e_k = zeros(1,length(x_bar)-1);
for i=2:length(x_bar)
    e_k = norm(x_bar(i,:) - x_bar(i-1,:));
    vec_e_k(1,i-1) = e_k;
end
figure;
idx_start = max(1, length(vec_e_k) - 200); 
plot(idx_start:length(vec_e_k), vec_e_k(idx_start:end), 'g', 'LineWidth', 2)
title('Study of  stagnation point','FontSize', 14)
xlabel('iteration', 'FontSize', 14)
ylabel('||x(k+1) - x(k)||', 'FontSize', 14)

