%% Wood function per dim = 10^3 %%
clc
clear all
close all
seed_value = min(343341, 343428);
rng(seed_value);
dim = 10^4;
disp("dimension:")
disp(dim)

% function definition
f2_wood = @(x) sum(100*(x(1:2:end-3).^2 - x(2:2:end-2)).^2 + (x(1:2:end-3) - 1).^2) + ...
    sum(90*(x(3:2:end-1).^2 - x(4:2:end)).^2 + (x(3:2:end-1) - 1).^2) + ...
    sum(10*(x(2:2:end-2) + x(4:2:end) - 2).^2 + (x(2:2:end-2) - x(4:2:end)).^2 / 10);

% gradient of the function:
function g = grad_wood_fun(x,h)
    n = length(x);
    g = sparse(n,1);

    g(1) = 400 * x(1)^3 - 400 * x(1) * x(2) + 2 * x(1) - 2;
    g(2) = -200 * x(1)^2 + (1101/5) * x(2) + 99/5 * x(4) - 40;

    g(n-1) = 360 * x(n-1)^3 - 360 * x(n-1) * x(n) + 2 * x(n-1) - 2;
    g(n) = (99/5) * x(n-2) - 180 * x(n-1)^2 + (1001/5) * x(n) - 40;
    
    % elementi 4:2:n-2
    g(4:2:n-2) = (99/5) * x(2:2:n-4) - 380 * x(3:2:n-3).^2 + (2102/5) * x(4:2:n-2) + (99/5) * x(6:2:n) - 80;
    
    % elementi 3:2:n-3
    g(3:2:n-3) = 760 * x(3:2:n-3).^3 + 4 * x(3:2:n-3) - 760 * x(3:2:n-3) .* x(4:2:n-2) - 4;
end
grad_wood = @(x,h) grad_wood_fun(x,h);

% Hessian matrix of the function:

% first implementation of the hessian matrix, when we realized it required
% too much time, we implemented the function above
% function H = hess_wood_fun_for_loops(x)
%     n = length(x);
%     H = sparse(n, n); 
% 
%     H(1,1) = 1200 * x(1)^2 - 400 * x(2) + 2;
%     H(2,2) = 1101/5;
%     H(n-1,n-1) = 1080 * x(n-1)^2 - 360 * x(n) + 2;
%     H(n,n) = 1001/5;
% 
%     H(1,2) = -400 * x(1);
%     H(n-1,n) = -360 * x(n-1);
% 
%     H(2,1) = H(1,2); 
%     H(n,n-1) = H(n-1,n);    
%     H(2,4) = 99/5;
%     H(n,n-2) = 99/5;
% 
%     for i = 3:2:n-3
%         H(i,i) = 2280 * x(i)^2 + 4 - 760 * x(i+1);
%         H(i,i+1) = -760 * x(i);
%     end
% 
%     for i = 4:2:n-2
%         H(i,i-2) = 99/5;
%         H(i,i-1) = -760 * x(i-1);
%         H(i,i) = 2102/5;
%         H(i,i+1) = 0;
%         H(i,i+2) = 99/5;
%     end
% end

function H = hess_wood_fun(x,h)
    n = length(x);

    % Pre-allocate vectors for diagonals
    main_diag = sparse(n,1);
    sup_diag1 = sparse(n,1); % sup-diagonal +1
    sup_diag2 = sparse(n,1); % sup-diagonal +2

    % --- Special elements --- %
    main_diag(1) = 1200 * x(1)^2 - 400 * x(2) + 2;
    sup_diag1(2) = -400 * x(1);
    main_diag(2) = 1101 / 5;
    sup_diag2(2) = 99 / 5;
    main_diag(n-1) = 1080 * x(n-1)^2 - 360 * x(n) + 2;
    sup_diag1(n) = -360 * x(n-1);
    main_diag(n) = 1001 / 5;

    % --- odd i = 3,5,...,n-3 --- %
    idx_odd = 3:2:n-3;
    main_diag(idx_odd) = 2280 * x(idx_odd).^2 + 4 - 760 * x(idx_odd+1);
    sup_diag1(idx_odd+1) = -760 * x(idx_odd);

    % --- even i = 4,6,...,n-2 --- %
    main_diag(idx_odd+1) = 2102 / 5;
    sup_diag2(idx_odd+3) = 99 / 5;

    sub_diag1 = [sup_diag1(2:n); 0];
    sub_diag2 = [sup_diag2(3:n); zeros(2,1)];

    % Final construction of the sparse hessian matrix
    H = spdiags([sub_diag2, sub_diag1, main_diag, sup_diag1, sup_diag2], ...
                [-2, -1, 0, 1, 2], n, n);
end
hess_wood = @(x,h) hess_wood_fun(x,h);

% Implementation of centered finite differences approximation for gradf and
% hessf:  
function g = gradf_finite_diff_wood_fun(x,h)
    n = length(x);
    g = sparse(n,1);

    g(1) = 400 * x(1)^3 - 400 * x(1) * x(2) + 2 * x(1) - 2 + 400 * h(1)^2 * x(1);
    g(2) = -200 * x(1)^2 + (1101/5) * x(2) + 99/5 * x(4) - 40;

    g(n-1) = 360 * x(n-1)^3 - 360 * x(n-1) * x(n) + 2 * x(n-1) - 2 + 180 * h(n-1) * x(n-1)^2;
    g(n) = (99/5) * x(n-2) - 180 * x(n-1)^2 + (1001/5) * x(n) - 40;
    
    % elements 4:2:n-2
    g(4:2:n-2) = (99/5) * x(2:2:n-4) - 380 * x(3:2:n-3).^2 + (2102/5) * x(4:2:n-2) + (99/5) * x(6:2:n) - 80;
    
    % elements 3:2:n-3
    g(3:2:n-3) = 760 * x(3:2:n-3).^3 + 4 * x(3:2:n-3) - 760 * x(3:2:n-3) .* x(4:2:n-2) - 4 + 760 * h(3:2:n-3) .* x(3:2:n-3);
end
gradf_finite_diff_wood = @(x,h) gradf_finite_diff_wood_fun(x,h);

function H = hess_finite_diff_wood_fun(x,h)
    n = length(x);

    % Prealloca vettori per le diagonali
    main_diag = sparse(n,1);
    sup_diag1 = sparse(n,1); % sopra diagonale a +1
    sup_diag2 = sparse(n,1); % sopra diagonale a +2

    % --- Elementi speciali ---
    main_diag(1) = 1200 * x(1)^2 - 400 * x(2) + 2 + 200 * h(1)^2;
    sup_diag1(2) = -400 * x(1) - 200 * h(1);
    main_diag(2) = 1101 / 5;
    sup_diag2(4) = 99 / 5; 
    main_diag(n-1) = 1080 * x(n-1)^2 - 360 * x(n) + 2 + 180 * h(n-1)^2;
    sup_diag1(n) = -360 * x(n-1) - 180 * h(n);
    main_diag(n) = 1001 / 5;

    % --- Loop sui dispari i = 3,5,...,n-3 ---
    idx_odd = 3:2:n-3;
    main_diag(idx_odd) = 2280 * x(idx_odd).^2 + 4 - 760 * x(idx_odd+1) + 380 * h(idx_odd).^2;
    sup_diag1(idx_odd+1) = -760 * x(idx_odd) - 380 * h(idx_odd);

    % --- Loop sui pari i = 4,6,...,n-2 ---
    main_diag(idx_odd+1) = 2102 / 5;
    sup_diag2(idx_odd+3) = 99 / 5;
    
    sub_diag1 = [sup_diag1(2:n); 0];
    sub_diag2 = [sup_diag2(3:n); zeros(2,1)];
    
    % Costruzione finale della matrice Hessiana sparsa
    H = spdiags([sub_diag2, sub_diag1, main_diag, sup_diag1, sup_diag2], ...
                [-2, -1, 0, 1, 2], n, n);
end
hess_finite_diff_wood = @(x,h) hess_finite_diff_wood_fun(x,h);

function H = hess_finite_diff_wood_h_variable(grad, x, h)
    n = length(x);
    H = sparse(n,n);

    % select gradient in positions: 1, 3, 5, 7, 9...
    jumps1 = repmat(2, 1, ceil(n/2));
    el1 = cumsum([1, jumps1]); % starts from position 1
    el1(el1 > n) = []; 
    x1 = x;
    x1(el1) = x(el1) + h(el1);
    x1_ = x;
    x1_(el1) = x(el1) - h(el1);

    % select gradient in positions: 2, 8, 14, 20...
    jumps2 = repmat(6, 1, ceil(n/2));
    el2 = cumsum([2, jumps2]); % starts from position 2
    el2(el2 > n) = []; 
    x2 = x;
    x2(el2) = x(el2) + h(el2);
    x2_ = x;
    x2_(el2) = x(el2) - h(el2);

    % select gradient in positions: 4, 10, 16, 22...
    el3 = cumsum([4, jumps2]); % starts from position 4
    el3(el3 > n) = []; 
    x3 = x;
    x3(el3) = x(el3) + h(el3);
    x3 = x;
    x3_(el3) = x(el3) - h(el3);

    % select gradient in positions: 6, 12, 18, 24...
    el4 = cumsum([6, jumps2]); % starts from position 6
    el4(el4 > n) = []; 
    x4 = x;
    x4(el4) = x(el4) + h(el4);
    x4_ = x;
    x4_(el4) = x(el4) - h(el4);

    col1 = [];
    h1 =  h(el1);
    h1_corr = [];
    for i = 1:length(el1)
        col1 = [col1, el1(i), el1(i)];
        h1_corr = [h1_corr, h1(i), h1(i)];
    end

    col2 = [];
    h2 =  h(el2);
    h2_corr = [];
    col2 = [col2, el2(1)*ones(1,5)];
    h2_corr = [h2_corr, h2(1)*ones(1,5)];
    for i = 2:length(el2)
        col2 = [col2, el2(i)*ones(1,6)];
        h2_corr = [h2_corr, h2(i)*ones(1,6)];
    end

    col3 = [];
    h3 =  h(el3);
    h3_corr = [];
    for i = 1:length(el3)
        col3 = [col3, el3(i)*ones(1,6)];
        h3_corr = [h3_corr, h3(i)*ones(1,6)];
    end

    col4 = [];
    h4 =  h(el4);
    h4_corr = [];
    for i = 1:length(el4)
        col4 = [col4, el4(i)*ones(1,6)];
        h4_corr = [h4_corr, h4(i)*ones(1,6)]; 
    end

    grad1 = grad(x1,h);
    grad2 = grad(x2,h);
    grad3 = grad(x3,h);
    grad4 = grad(x4,h);
    grad1_ = grad(x1_,h);
    grad2_ = grad(x2_,h);
    grad3_ = grad(x3_,h);
    grad4_ = grad(x4_,h);

    % H(1,col1(1)) = (grad1(1) - gradx(1))./(h1_corr(1));
    % H(1,col2(1)) = (grad2(1) - gradx(1))./(h2_corr(1));
    H(1,col1(1)) = (grad1(1) - grad1_(1))./(2*h1_corr(1));
    H(1,col2(1)) = (grad2(1) - grad2_(1))./(2*h2_corr(1));
    for i=2:n
        % H(i,col1(i)) = (grad1(i) - gradx(i))./(h1_corr(i));
        % H(i,col2(i)) = (grad2(i) - gradx(i))./(h2_corr(i));
        % H(i,col3(i-1)) = (grad3(i) - gradx(i))./(h3_corr(i-1));
        H(i,col1(i)) = (grad1(i) - grad1_(i))./(2*h1_corr(i));
        H(i,col2(i)) = (grad2(i) - grad2_(i))./(2*h2_corr(i));
        H(i,col3(i-1)) = (grad3(i) - grad3_(i))./(2*h3_corr(i-1));
        if i >= 4 && i <= n-2
            % H(i,col4(i-3)) = (grad4(i) - gradx(i))./(h4_corr(i-3));
            H(i,col4(i-3)) = (grad4(i) - grad4_(i))./(2*h4_corr(i-3));
            % set to 0 components which are 0
            H(i-1,i-2) = 0;
        end
    end

    H = 0.5 * (H + H');
end
hess_finite_diff_wood_h_var = @(grad, x,h) hess_finite_diff_wood_h_variable(grad, x, h);

n = 1:dim; 
x2_wood = zeros(dim,1); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;  
x2_opt = ones(dim,1);
l_bound = x2_wood - 1;
u_bound = x2_wood + 1; 
M_ten_initial_points = l_bound + (u_bound - l_bound) .* rand(dim, 10);
starting_points = [x2_wood, M_ten_initial_points];

print_function(x2_opt, f2_wood)

% obs: rho = 0.9 and c1 = 10^(-4), btmax = 50 -> changes too slow
% obs: rho = 0.9 and c1 = 10^(-4), btmax = 100 -> changes too slow
% rho = 0.1 % too low
rho = 0.9; 
tolgrad = 10^(-2);
c1 = 10^(-3); 
btmax = 90;
kmax = 700;

flag_h = 0;
h = 0;

diff_wood = zeros(11,1);
time_wood = zeros(11,1);
iter_wood = zeros(11,1);
rate_wood = zeros(1,11);
failure_wood = zeros(1,11);
for index=1:11
    tic;
    [xk, fk, gradfk_norm, k, ...
        xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f2_wood, grad_wood, hess_wood, kmax, ...
                                                            tolgrad, c1, btmax, dim, rho, h, flag_h);
    time_wood(index) = toc;
    iter_wood(index) = k;
    diff_wood(index) = norm(xk - x2_opt,2);
    rate_wood(index) = compute_exp_rate_conv_multi(xseqk);
    failure_wood(index) = failure;        
    disp(gradfk_norm)
end
final_rate_exact_deriv_wood = mean(rate_wood, 'omitnan');

% Finite differences

% parameter used for finite differences
h_vec = 10.^[-2, -4, -6, -8, -10, -12];

% Testing different values of h, fixed and varying from point to point 
flag_h2 = 1; % h must be adapted to each component of x

num_point = 11;

time1_finite_diff_wood = zeros(length(h_vec),num_point);
rate_h_wood = zeros(1, num_point);
final_rate_wood = zeros(1, length(h_vec));
norm_err_conv_wood = zeros(num_point, length(h_vec));
iters_wood = zeros(length(h_vec), num_point);
failure1_wood = zeros(length(h_vec),num_point);

time2_finite_diff_wood = zeros(length(h_vec),num_point);
rate_h2_wood = zeros(1, num_point);
final_rate2_wood = zeros(1,length(h_vec));
norm_err_conv2_wood = zeros(num_point, length(h_vec));
iters2_wood = zeros(length(h_vec), num_point);
failure2_wood = zeros(length(h_vec),num_point);

final_point_temporary_wood = zeros(dim,num_point);
final_point_wood = cell(1,num_point);
final_point_temporary2_wood = zeros(dim,num_point);
final_point2_wood = cell(1,num_point);

% obs dim=10^3: for all initial point when h = 10^-2, 10^-3 it is displayed the warning:
% Warning: Couldn't find an optimal alpha for bactracking. 
% obs dim=10^4: for all initial point when h = 10^-2, 10^-3 it is displayed the warning:
% Warning: Couldn't find an optimal alpha for bactracking. 
% obs dim=10^5: for h=10^-2 and initial points number 1,2,4,5 it is displayed:
% Warning: Couldn't find an optimal alpha for bactracking. 
% obs: increasing the parameter btmax permits to diminuish the warning:
% Warning: Couldn't find an optimal alpha for bactracking. 
for i=1:length(h_vec)
    disp(i)
    for index=1:num_point
         disp(index)
        % tic;
        % [xk, fk, gradfk_norm, k, ...
        %     xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f2_wood, ...
        %                                                     gradf_finite_diff_wood, hess_finite_diff_wood, kmax,...
        %                                                     tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h);
        % time1_finite_diff_wood(i, index) = toc;
        % iters_wood(i,index) = k;
        % norm_err_conv_wood(index,i) = norm(x2_opt - xk,2);
        % rate_h_wood(index) = compute_exp_rate_conv_multi(xseqk);
        % final_point_temporary_wood(:,index) = xk;
        % failure1_wood(i,index) = failure;
        % disp(rate_h_wood(index))
        % disp(gradfk_norm)

        tic;
        [xk2, fk2, gradfk_norm2, k2, ...
            xseqk2, btseq2, failure2] = modified_newton_method(starting_points(:,index), f2_wood,...
                                                  grad_wood, hess_finite_diff_wood, kmax,...
                                                  tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h2);
        time2_finite_diff_wood(i, index) = toc;
        iters2_wood(i, index) = k2;
        norm_err_conv2_wood(i,index) = norm(x2_opt - xk2,2);
        rate_h2_wood(index) = compute_exp_rate_conv_multi(xseqk2);
        final_point_temporary2_wood(:,index) = xk2;
        failure2_wood(i,index) = failure2;
    end
    final_rate_wood(i) = mean(rate_h_wood, 'omitnan');
    final_rate2_wood(i) = mean(rate_h2_wood, 'omitnan');
    final_point_wood{i} = final_point_temporary_wood;
    final_point2_wood{i} = final_point_temporary2_wood;
end

% filename = 'output.mat';
% 
% %base_vars = {'time_wood', 'diff_wood', 'iter_wood', 'final_rate_exact_deriv_wood', 'failure_wood'};
% base_vars = {'time1_finite_diff_wood', 'norm_err_conv_wood', 'iters_wood', 'final_rate_wood', 'failure1_wood'};
% %base_vars = {'time2_finite_diff_wood', 'norm_err_conv2_wood', 'iters2_wood', 'final_rate2_wood', 'failure2_wood'};
% 
% for i = 1:length(base_vars)
%     base_name = base_vars{i};
%     new_name = sprintf('%s_%d', base_name, dim);  % es. time_pow_1000
% 
%     val = eval(base_name);
% 
%     eval([new_name ' = val;']);
% 
%     if exist(filename, 'file')
%         save(filename, new_name, '-append');
%     else
%         save(filename, new_name);
%     end
% end
% 
