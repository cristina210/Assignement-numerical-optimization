%% Chained Rosenbrook for dim = 10^3 %%
clc
clear all
close all

% set the seed to obtain the same results while running the program
seed_value = min(343341, 343428);
rng(seed_value);
dim = 10^4;
disp("dimension:")
disp(dim)

% function definition
f1_ros = @(x) sum(100*(x(1:end-1).^2 - x(2:end)).^2 + (x(1:end-1) - 1).^2);

% gradient of the function
grad_f1_ros = @(x,h) sparse([
    2 * x(1) + 400 * x(1)^3 - 400 * x(1)*x(2) - 2;
    -200 * x(1:end-2).^2 + 202 * x(2:end-1) + 400 * x(2:end-1).^3 - 400 * x(2:end-1).*x(3:end) - 2;
    -200 * x(end - 1)^2 + 200 * x(end) 
]);  

% Hessian matrix of the function
hess_f1_ros = @(x,h) spdiags([
    sparse([-400*x(1:end-1); 0]), ...
    sparse([1200 * x(1)^2 - 400 * x(2) + 2; ...
            202 + 1200 * x(2:end - 1).^2 - 400 * x(3:end); ...
            200]), ...
    sparse([0; -400*x(1:end-1)])
], [-1, 0, 1], length(x), length(x));


% Implementation of centered finite differences approximation for gradf and
% hessf: 
gradf_finite_diff_ros = @(x,h) [400*x(1)^3 + 400*h(1)^2*x(1) - 400*x(1)*x(2) + 2*x(1) - 2; ...
    -200*x(1:end-2).^2 + 202*x(2:end-1) + 400*x(2:end-1).^3 + 400*h(2:end-1).^2.*x(2:end-1) - 400*x(2:end-1).*x(3:end) - 2; ...
    -200*x(end-1)^2 + 200*x(end)];


hess_finite_diff_ros = @(x, h) spdiags([
    [-400 * x(1:end-1) - 200 * h; 0],...
    [1200 * x(1)^2 + 200 * h^2 - 400 * x(2) + 2; ...
    1200 * x(2:end-1).^2 - 400 * x(3:end) + 200 * h^2 + 202; ...
    200],...
    [0; -400 * x(1:end-1) - 200 * h]
    ], [-1, 0, 1 ], length(x), length(x));

% Function computed exploiting the least number of function evaluations
% function H = hess_finite_diff_ros_h_variable(grad, x, h)
%     n = length(x);
%     H = sparse(n,n);
%     gradx = grad(x);
% 
%     % select gradient in positions: 1, 4, 7...
%     jumps1 = repmat(3, 1, ceil(n/2));
%     el1 = cumsum([1, jumps1]); % starts from position 1
%     el1(el1 > n) = []; 
%     x1 = x;
%     x1(el1) = x(el1) + h(el1);
%     x1_ = x;
%     x1_(el1) = x(el1) - h(el1);
% 
%     % select gradient in positions: 2, 5, 8...
%     el2 = cumsum([2, jumps1]); % starts from position 2
%     el2(el2 > n) = []; 
%     x2 = x;
%     x2(el2) = x(el2) + h(el2);
%     x2_ = x;
%     x2_(el2) = x(el2) - h(el2);
% 
%     % select gradient in positions: 3, 6, 9...
%     el3 = cumsum([3, jumps1]); % starts from position 3
%     el3(el3 > n) = []; 
%     x3 = x;
%     x3(el3) = x(el3) + h(el3);
%     x3_ = x;
%     x3_(el3) = x(el3) - h(el3);
% 
%     col1 = [];
%     h1 =  h(el1);
%     h1_corr = [];
% 
%     col1 = [el1(1), el1(1), repelem(el1(2:end), 3)];
%     h1_corr = [h1(1), h1(1), repelem(h1(2:end), 3)];
% 
%     % for i = 2:length(el1)
%     %     if i == 1
%     %         col1 = [col1, el1(1), el1(1)];
%     %         h1_corr = [h1_corr, h1(1), h1(1)];
%     %     else
%     %     col1 = [col1, el1(i),  el1(i),  el1(i)];
%     %     h1_corr = [h1_corr, h1(i), h1(i), h1(i)];
%     %     end
%     % end
% 
%     col2 = [];
%     h2 =  h(el2);
%     h2_corr = [];
%     % for i = 1:length(el2)
%     %     col2 = [col2, el2(i), el2(i), el2(i)];
%     %     h2_corr = [h2_corr, h2(i), h2(i), h2(i)];
%     % end
%     col2 = repelem(el2, 3);
%     h2_corr = repelem(h2, 3);
% 
%     col3 = [];
%     h3 =  h(el3);
%     h3_corr = [];
%     % for i = 1:length(el3)
%     %     col3 = [col3, el3(i), el3(i), el3(i)];
%     %     h3_corr = [h3_corr, h3(i), h3(i), h3(i)];
%     % end
%     col3 = repelem(el3, 3);
%     h3_corr = repelem(h3, 3);
% 
% 
%     grad1 = grad(x1);
%     grad1_ = grad(x1_);
%     grad2 = grad(x2);
%     grad2_ = grad(x2_);
%     grad3 = grad(x3);
%     grad3_ = grad(x3_);
% 
%     col2 = [col2, 0];
%     col3 = [0, col3];
% 
%     H(1,col1(1)) = (grad1(1) - gradx(1))./(2*h1_corr(1));
%     H(1,col2(1)) = (grad2(1) - gradx(1))./(2*h2_corr(1));
%     for i=2:length(x)-1
%         H(i,col1(i)) = (grad1(i) - gradx(i))./(2*h1_corr(i));
%         H(i,col2(i)) = (grad2(i) - gradx(i))./(2*h2_corr(i));
%         H(i,col3(i)) = (grad3(i) - gradx(i))./(2*h3_corr(i));
%     end
%     H(end,col1(end)) = (grad1(end) - gradx(end))./(2*h1_corr(end));
%     H(end,col3(end)) = (grad3(end) - gradx(end))./(2*h3_corr(end));
% 
%     H = 0.5 * (H + H');
% end

hess_finite_diff_ros_h_var = @(x, h) spdiags([
    [-400 * x(1:end-1) - 200 * h(1:end-1); 0],...
    [1200 * x(1)^2 + 200 * h(1)^2 - 400 * x(2) + 2; ...
    1200 * x(2:end-1).^2 - 400 * x(3:end) + 200 * h(2:end-1).^2 + 202; ...
    200],...
    [0; -400*x(1:end-1) - 200 * h(1:end-1)]
    ], [-1, 0, 1 ], length(x), length(x));

x1_rosenbrock = repmat([-1.2; 1.0], dim/2, 1); % because dim is even
x1_opt = ones(dim,1); 
l_bound = x1_rosenbrock - 1;
u_bound = x1_rosenbrock + 1; 
M_ten_initial_points = l_bound + (u_bound - l_bound) .* rand(dim, 10);
starting_points = [x1_rosenbrock, M_ten_initial_points];

print_function(x1_opt, f1_ros)

rho = 0.9;
tolgrad = 10^(-2);
c1 = 10^(-4); 
btmax = 75;
% kmax = 1500;  % for dim = 10^3
% kmax = 15000; % for dim = 10^4
kmax = 150000;  % for dim = 10^5

flag_h = 0;
h = 0;

diff_ros = zeros(11,1);
time_ros = zeros(11, 1);
iter_ros = zeros(11, 1);
rate_ros = zeros(1,11); 
failure_ros = zeros(1,11);
for index=1:11
    tic;
    [xk, fk, gradfk_norm, k, ...
        xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f1_ros, grad_f1_ros, hess_f1_ros, kmax, ...
                                                            tolgrad, c1, btmax, dim, rho, h, flag_h);
    time_ros(index) = toc;
    iter_ros(index) = k;
    diff_ros(index) = norm(xk - x1_opt,2);
    rate_ros(index) = compute_exp_rate_conv_multi(xseqk);
    failure_ros(index) = failure;
    disp(gradfk_norm)
end
final_rate_exact_deriv_ros = mean(rate_ros, 'omitnan');

% Finite differences
% parameter used for finite differences
h_vec = 10.^[-2, -4, -6, -8, -10, -12];

% Testing different values of h, fixed and varying from point to point 
flag_h2 = 1; % Hessian matrix computed with finite differences and h variable

num_point = 11;

time1_finite_diff_ros = zeros(length(h_vec),num_point);
rate_h_ros = zeros(1, num_point);
final_rate_ros = zeros(1, length(h_vec));
norm_err_conv_ros = zeros(num_point, length(h_vec));
iters_ros = zeros(num_point, length(h_vec));
failure1_ros = zeros(length(h_vec),num_point);

time2_finite_diff_ros = zeros(length(h_vec),num_point);
rate_h2_ros = zeros(1, num_point);
final_rate2_ros = zeros(1,length(h_vec));
norm_err_conv2_ros = zeros(num_point, length(h_vec));
iters2_ros = zeros(num_point, length(h_vec));
failure2_ros = zeros(length(h_vec),num_point);

final_point_temporary_ros = zeros(dim,11);
final_point_ros = cell(1,11);
final_point_temporary2_ros = zeros(dim,11);
final_point2_ros = cell(1,11);

for i=1:length(h_vec)
    disp(i)
    for index=1:num_point
        %disp(index)
        % tic;
        % [xk, fk, gradfk_norm, k, ...
        %     xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f1_ros, ...
        %                                                     gradf_finite_diff_ros, hess_finite_diff_ros_h_var, kmax,...
        %                                                     tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h);
        % time1_finite_diff_ros(i, index) = toc;
        % iters_ros(index, i) = k;
        % norm_err_conv_ros(index,i) = norm(x1_opt - xk,2);
        % rate_h_ros(index) = compute_exp_rate_conv_multi(xseqk);
        % final_point_temporary_ros(:,index) = xk;
        % failure1_ros(i, index) = failure;

        disp(index)
        tic;
        [xk2, fk2, gradfk_norm2, k2, ...
            xseqk2, btseq2, failure2] = modified_newton_method(starting_points(:,index), f1_ros,... 
                                                               grad_f1_ros, hess_finite_diff_ros_h_var, kmax,...
                                                               tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h2);
        time2_finite_diff_ros(i, index) = toc;
        iters2_ros(index, i) = k2;
        norm_err_conv2_ros(index,i) = norm(x1_opt - xk2,2);
        rate_h2_ros(index) = compute_exp_rate_conv_multi(xseqk2);
        final_point_temporary2_ros(:,index) = xk2;
        failure2_ros(i, index) = failure2;
    end
    final_rate_ros(i) = mean(rate_h_ros, 'omitnan');
    final_rate2_ros(i) = mean(rate_h2_ros, 'omitnan');
    final_point_ros{i} = final_point_temporary_ros;
    final_point2_ros{i} = final_point_temporary2_ros;
end


filename = 'output.mat';
% base_vars = {'time_ros', 'diff_ros', 'iter_ros', 'final_rate_exact_deriv_ros', 'failure_ros'};
% base_vars = {'time1_finite_diff_ros', 'norm_err_conv_ros', 'iters_ros', 'final_rate_ros', 'failure1_ros'};
base_vars = {'time2_finite_diff_ros', 'norm_err_conv2_ros', 'iters2_ros', 'final_rate2_ros', 'failure2_ros'};

for i = 1:length(base_vars)
    base_name = base_vars{i};
    new_name = sprintf('%s_%d', base_name, dim);  

    val = eval(base_name);

    eval([new_name ' = val;']);

    if exist(filename, 'file')
        save(filename, new_name, '-append');
    else
        save(filename, new_name);
    end
end

