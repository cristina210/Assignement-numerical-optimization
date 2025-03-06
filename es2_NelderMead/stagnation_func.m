% Function: stagnation_func
% This function calculates and visualizes the norm of the differences
% between consecutive solutions. It is intended to 
% study the presence of a "stagnation point" by plotting the changes in 
% successive iterations.
%
% Inputs:
% - x: A matrix where each row represents an iteration point, and 
%   columns represent value of the solution's coordinates.
%
% Outputs:
% - vec_e_k: A row vector containing norms of the differences 
%   between consecutive solutions.

function [vec_e_k] = stagnation_func(x)
vec_e_k = zeros(1,length(x)-1);
for i=2:length(x)
    e_k = norm(x(i,:) - x(i-1,:));
    vec_e_k(1,i-1) = e_k;
end
figure;
plot(1:length(vec_e_k), vec_e_k, 'g', 'LineWidth', 2)
title('Study of  stagnation point','FontSize', 16)
xlabel('iteration', 'FontSize', 14)
ylabel('||x(k+1) - x(k)||', 'FontSize', 14)

