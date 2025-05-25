function [xk, fk, gradfk_norm, k, ...
    xseq, btseq, failure] = modified_newton_method(x0, f, gradf, hessf, kmax, tolgrad, c1, btmax, n, rho, h0, flag_h)

% function [xk, fk, gradfk_norm, k, xseq, btseq, failure] = modified_newton_method(x0, f, gradf, hessf, kmax, tolgrad, c1, btmax, n, rho, h0, flag_h)
% Function that performs the Newton optimization method when the hessian
% matrix is not positive definite, using bactracking strategy for the
% step-length selection.
% 
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t. the norm of the
% gradient;
% c1 = ﻿the factor of the Armijo condition that must be a scalar in (0,1);
% btmax = ﻿maximum number of steps for updating alpha during the 
% backtracking strategy;
% n = dimension of x0;
% rho = ﻿fixed factor, lesser than 1, used for reducing alpha;
% h0 = step-length through which we compute the expansion x = x_bar + h used
% to compute finite difference;
% flag_h = flag used to indicate which h to use, if flag_h = 1, need to
% adapt h at each value, if flag_h = 0, h is costant, even in the case with
% exact derivatives where h = 0;
% 
% OUTPUTS:
% xk = the last x computed by the function with preconditioning;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk);
% k = index of the last iteration performed;
% xseq = 4-by-k matrix where the columns are the elements xk of the 
% sequence;
% btseq = 1-by-k vector where elements are the number of backtracking
% iterations at each optimization step with xk;
% failure = flag which indicates if there is stagnation during
% backtracking.


% Pre-allocate
xk = x0;
fk = f(xk);
xseq = zeros(n, 4);
btseq = zeros(1, kmax);
beta = 10^-3;
failure = false;
h = ones(n,1) * h0;

if flag_h == 1 
    h = h0*abs(xk);
end

gradfk = gradf(xk,h);
gradfk_norm = norm(gradfk);

farmijo = @(fk, alpha, c1_gradfk_pk) fk + alpha * c1_gradfk_pk;
k = 0;

% best_values = zeros(kmax,1);
% best_values(1) = fk;
% best_gradf = zeros(kmax,1);
% best_gradf(1) = norm(gradfk);

time_limit = 1000;
start_time = tic;
while k < kmax && gradfk_norm > tolgrad

    if toc(start_time) > time_limit
        failure = true;
        warning('The method stopped because it reached time_limit')
        break;
    end

    if flag_h == 1 
        h = h0*abs(xk);
    end
    hessfk = hessf(xk,h);

    % solve with Cholesky correction
    L = choleski_added_multiple_identity(hessfk, beta, n);
    y = L \ (-gradfk);
    pk = L' \ y;

    % Backtracking
    alpha = 1;
    xnew = xk + alpha * pk;
    fnew = f(xnew);
    c1_gradfk_pk = c1 * gradfk' * pk;
    bt = 0;

    while bt < btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        alpha = rho * alpha;
       
        xnew = xk + alpha * pk;

        if norm(xnew-xk,2) < 10^-12
            failure = true;
            warning('Stagnation in backtracking')
            break
        end
        fnew = f(xnew);
        bt = bt + 1;
    end

    if bt == btmax && fnew > farmijo(fk, alpha, c1_gradfk_pk)
        warning('Backtracking failed to find a good step.');
        failure = true;
        btseq(k+1, 1) = bt; 
        xk = xnew;
        k = k+1;
        break;
    end

    if failure == true
        warning("Stagnation in bactracking: change btmax or rho parameters.")
        break
    end

    % Update
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk,h);
    gradfk_norm = norm(gradfk);

    k = k + 1;
    xseq(:, mod(k, 4) + 1) = xk;
    btseq(k) = bt;


    % best_values(k) = fk;
    % best_gradf(k) = norm(gradfk);
    % if mod(k, 10) == 0
    %     figure(2);
    %     plot(best_values(6:k), '-o', 'MarkerSize', 6, 'LineWidth', 2);
    %     xlabel('Iterations', 'FontSize', 14);
    %     ylabel('Best Evaluation', 'FontSize', 14);
    %     title('Progress minimum value Modified Newton Method', 'FontSize', 16);
    %     set(gca, 'FontSize', 12); 
    %     drawnow;
    % 
    %     figure(3);
    %     plot(best_gradf(5:k), '-o', 'MarkerSize', 6, 'LineWidth', 3);
    %     xlabel('Iterations', 'FontSize', 14);
    %     ylabel('Best Evaluation', 'FontSize', 14);
    %     title('Progress gradient value Modified Newton Method', 'FontSize', 16);
    %     set(gca, 'FontSize', 12); 
    %     drawnow;
    % end

end

if (k == kmax || failure) && gradfk_norm > tolgrad
    failure = true;
end

btseq = btseq(1:k);

% Add final block
% if mod(k, 10) ~= 0 && k > 0
%     figure(1);
%     %hold off;
%     plot(best_values(1:k), '-o', 'MarkerSize', 6, 'LineWidth', 2);
%     xlabel('Iterations', 'FontSize', 14);
%     ylabel('Best Evaluation', 'FontSize', 14);
%     title('Progress minimum value Modified Newton Method', 'FontSize', 16);
%     set(gca, 'FontSize', 12);
%     drawnow;
% end


end