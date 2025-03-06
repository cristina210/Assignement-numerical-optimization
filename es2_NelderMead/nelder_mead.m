% Function: nelder_mead
% Implements the Nelder-Mead method to minimize a given function
%
% Inputs:
% - f: function  to minimize.
% - simplex: initial simplex (rows = points, cols = coordinate).
% - kmax: maximum number of iterations.
% - rho, chi, gamma, sigma:  coefficients for reflection, expansion, 
%   contraction, and shrinkage.
% - n: dimension
% - tol1: Tolerance for the stopping criterion based on the size of the 
%   simplex, indicating convergence to a single point.
% - tol2: Tolerance for the stopping criterion based on the function value, 
%   triggered when the simplex reaches a stationary region for f.
% Outputs:
% - k: number of iterations performed.
% - simplex: final simplex.
% - x_vec: vector with simplex best point during iterations.
% - flag: Status flag indicating the output of the algorithm
%   (0: no issue, 1: max iterations reached, 2: stopping criterion based on tol1 is
%    met, 3: stopping criterion based on tol2 is met).
% - type_of_move: is a matrix that tracks the type of transformation applied to 
%   the simplex at each iteration of the Nelder-Mead algorithm. Each row contains a triplet. 
% - size_vec: vector with the max distance simplex's vertice - barycenter
%   through iterations


function [k, simplex, x_vec, flag, type_of_move, size_vec] = nelder_mead(f, simplex, kmax, rho, chi, gamma, n, sigma, tol1, tol2)
flag = 0;
tol3 = 10^(-8); 
type_of_move = [];
%[vol, flagg] = volume_simplex(simplex);

% Function handle for updating certain quantities
upDate_quantities = @(k, simplex, x_bar, f1, fend) deal( ...
    k+1, ...
    max(vecnorm(simplex - x_bar, Inf, 2)), ...
    abs(fend - f1));

% initialization
history = [];
 
f_val = zeros(n+1, n+1);
for i=1:n+1
    f_val(i,:) = [f(simplex(i,:)), simplex(i,:)];
end
% f_val is a matrix of size (n+1, n+1) that contains the function values 
% and the coordinates of the points in the simplex during the execution 
% of the Nelder-Mead algorithm. The first column of f_val contains the function 
% values f computed at the points of the simplex. 
% Each row represents a point in the simplex, 
% so f_val(i,1) is the function value evaluated at the i-th point.

x_bar = mean(simplex(1:n, :));

x_vec = simplex(1,:);
% x_vec contains best point of the simplex through iterations

distance_bar = max(vecnorm(simplex - x_bar, Inf, 2));
% distance_bar represents the maximum distance between the barycenter (x_bar) 
% and the points of the simplex.
% distance_bar is used as a convergence criterion: when it becomes smaller than 
% a predefined tolerance (tol1) the method stops.

size_vec = [distance_bar];
% size_vec contains the maximum distance between a vertice of the symplex
% anche the barycenter rapresenting the dimension of the simplex

delta_f = 1;
% delta_f represents the absolute difference between the function values 
% at the first and last points of the simplex.
% delta_f is used as a convergence criterion: when it becomes smaller than 
% a predefined tolerance (tol2) the method stops.

k = 0;


while k < kmax && distance_bar > tol1 && delta_f > tol2
    if k ~= 1
        x_vec = [x_vec; simplex(1,:)];
    end
    f_val = sortrows(f_val);
    % f_val is sorted based on the function values to ensure 
    % that the best point is always at the top of the simplex.

    simplex = f_val(:,2:end); 
    
    x_bar = mean(simplex(1:n, :));

    % Save the volume of simplex
    %[volume, flagg] = volume_simplex(simplex);    

    % reflection
    x_r = x_bar + rho*(x_bar - simplex(n+1,:));
    f1 = f_val(1,1);
    fr = f(x_r); 
    fn = f_val(n, 1);
    fend = f_val(n+1, 1);
    if (f1 < fr || abs(f1-fr) <= tol3) && ( fr < fn || abs(fr-fn) < tol3)
        % enough reduction but not exceptional
        simplex(end,:) = x_r;
        f_val(end,:) = [fr, x_r];
        [k, distance_bar, delta_f] = upDate_quantities(k, simplex, x_bar, f1, fend);
        type_of_move = [type_of_move; [1,0,1]];  % colour magenta
        size_vec = [size_vec; distance_bar];
        continue 
    elseif fr < f1 
        % expansion: exceptional reduction
        x_e = x_bar + chi*(x_r - x_bar);
        fe = f(x_e);
        history = [history; "exp"];
        if fe < fr
            simplex(end, :) = x_e; 
            f_val(end,:) = [fe, x_e];
            [k, distance_bar, delta_f] = upDate_quantities(k, simplex, x_bar, f1, fend);
            type_of_move = [type_of_move; [0,1,0]];  % colour green 
            size_vec = [size_vec; distance_bar];
            continue
        else
            simplex(end,:) = x_r;            
            f_val(end,:) = [fr, x_r];
            [k, distance_bar, delta_f] = upDate_quantities(k, simplex, x_bar, f1, fend);
            type_of_move = [type_of_move; [0,1,0]];  % colour green 
            size_vec = [size_vec; distance_bar];
            continue
        end
    elseif (fr > fn || abs(fn-fr) <= tol3)
        % contraction: absente or very little reduction
        if fend < fr
            x_c = x_bar - gamma*(x_bar - simplex(end,:));
        else
            x_c = x_bar - gamma*(x_bar - x_r);
        end
        fc = f(x_c);
        history = [history; "cont"];
        if fc < fend 
            % new point is better than the worst point in the simplex
            simplex(end, :) = x_c;
            f_val(end,:) = [fc, x_c];
            [k, distance_bar, delta_f] = upDate_quantities(k, simplex, x_bar, f1, fend);
            type_of_move = [type_of_move; 	[1,1,0]];  % colour yellow
            size_vec = [size_vec; distance_bar];
            continue
        else
            % worst situation
            % shrinking
            for i=2:n+1
                simplex(i,:) = simplex(1,:) + sigma*(simplex(i,:) - simplex(1,:));
                f_val(i,:) = [f(simplex(i,:)), simplex(i,:)];
            end
            history = [history; "shrink"];
            [k, distance_bar, delta_f] = upDate_quantities(k, simplex, x_bar, f1, fend);
            type_of_move = [type_of_move; [1,0,0]];  % colour red
            size_vec = [size_vec; distance_bar];
            continue
        end
    end
end
count_shrinks = sum(history == "shrink");
count_contr = sum(history == "cont");
count_exp = sum(history == "exp");
disp("How many shrink:")
disp(count_shrinks)
disp("How many expansion:")
disp(count_exp)
disp("How many contraction:")
disp(count_contr)
x_vec = [x_vec; simplex(1,:)];
x_vec = x_vec(2:size(x_vec,1),:);

if k == 1
   warning('The initial simplex provided as input is not adequate to explore the surrounding space: the function stops at the first iteration')
end
if k == kmax
   warning('Maximum iteration limit reached')
   flag = 3;  % no convergence
end
if distance_bar <= tol1 
   disp('The simplex is converging to a point')
   flag = 1; % convergence of type 1
end
if delta_f <= tol2 
   disp('The simplex has reached a stationary region for f')
   flag = 2; % convergence of type 2
end

end