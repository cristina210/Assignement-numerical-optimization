% Function: NelderMead_simplex
% Constructs an initial simplex for the Nelder-Mead method 
% in a specified dimension. Verifies that the simplex is non-degenerate 
% by checking the linear independence of its vertices.
%
% Inputs:
% - dim: dimension
% - initial_point: coordinates of a point for the simplex.
%
% Outputs:
% - simplex_initial: Matrix containing the vertices of the initial simplex.
% - flag: Status flag indicating whether the simplex is valid.

function [simplex_initial, flag] = NelderMead_simplex(dim, initial_point)
flag = 0;

% Create initial simplex
simplex_initial = zeros(dim+1, dim); 
for i = 1:(dim+1)
    for j = 1:(dim)
        if i >= 2 && i == j + 1
            simplex_initial(i,j) = 1;
        else
            simplex_initial(1,:) = initial_point;
        end
    end
end

% Check if it's degenerate
vettore_differenze = zeros(dim,dim);
for i = 2:dim+1
    vettore_differenze(i-1,:) = simplex_initial(1,:) - simplex_initial(i,:);
end
if rank(vettore_differenze) ~= dim
    disp("Initial symplex is invalid")
    flag = 1;
end
end

