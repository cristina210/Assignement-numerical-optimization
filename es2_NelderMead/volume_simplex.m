function [volume, flag] = volume_simplex(simplex)
dim = size(simplex, 2);
flag = 0; 

mat_difference = zeros(dim, dim);
for i = 2:dim+1
    mat_difference (i-1, :) = simplex(i, :) - simplex(1, :);
end

% Check if the simplex is valid (non-degenerate)
if rank(mat_difference ) ~= dim
    disp("Initial simplex is invalid")
    flag = 1;   % simplex is invalid
    V = 0; 
    return;
end

% Compute the volume of the simplex
volume = abs(det(mat_difference )) / factorial(dim);
end

