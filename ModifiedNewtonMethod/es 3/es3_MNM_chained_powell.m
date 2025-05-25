%% Powell function per dim = 10^5 %%
clc
clear all
close all
seed_value = min(343341, 343428);
rng(seed_value);
dim = 10^5;
disp("dimension:")
disp(dim)

% function definition
f3_powell = @(x) sum(arrayfun(@(j) ...
        (x(2*j-1) + 10*x(2*j))^2 + ...
        5*(x(2*j+1) - x(2*j+2))^2 + ...
        (x(2*j) - 2*x(2*j+1))^4 + ...
        10*(x(2*j-1) - x(2*j+2))^4, ...
        1:(length(x)-2)/2));

% gradient of the function
function g = grad_powell_fun(x,h)
    n = length(x);
    g = sparse(n,1); 

    g(1) = 40 * x(1)^3 + 2 * x(1) + 20 * x(2) - 120 * x(1)^2 * x(4) + 120 * x(1) * x(4)^2 - 40 * x(4)^3; 
    g(2) = 20 * x(1) + 200 * x(2) + 4 * x(2)^3 - 24 * x(2)^2 * x(3) + 48 * x(2) * x(3)^2 - 32 * x(3)^3; 
 
    idx_odd = 3:2:n-3;

    g(idx_odd) = 12 * x(idx_odd) + 10 * x(idx_odd-1) - 8 * x(idx_odd-1).^3 + 48 * x(idx_odd-1).^2 .* x(idx_odd) ...
        - 96 * x(idx_odd-1) .* x(idx_odd).^2 + 104 * x(idx_odd).^3 - 120 * x(idx_odd).^2 .* x(idx_odd+3) + 120 * x(idx_odd) .* x(idx_odd+3).^2 ...
        - 40 * x(idx_odd+3).^3; 

    g(idx_odd+1) = 10 * x(idx_odd) + 210 * x(idx_odd+1) - 40 * x(idx_odd-2).^3 + ...
        120 * x(idx_odd-2).^2 .* x(idx_odd+1) - 120 * x(idx_odd-2) .* x(idx_odd+1).^2 + ...
        44 * x(idx_odd+1).^3 - 24 * x(idx_odd+1).^2 .* x(idx_odd+2) + 48 * x(idx_odd+1) .* x(idx_odd+2).^2 - 32 * x(idx_odd+2).^3; 
    

    g(n-1) = 10 * x(n-1) - 10 * x(n) - 8 * x(n-2)^3 + 48 * x(n-2)^2 * x(n-1) - 96 * x(n-2) * x(n-1)^2 + 64 * x(n-1)^3; 
    g(n) = - 10 * x(n-1) + 10 * x(n) - 40 * x(n-3)^3 + 120 * x(n-3)^2 * x(n) - 120 * x(n-3) * x(n)^2 + 40 * x(n)^3;
end
grad_powell = @(x,h) grad_powell_fun(x,h);

% Hessian matrix of the function:
function H = hess_powell_fun(x,h)
    n = length(x);
     
    % Pre-allocate vectors for diagonals
    main_diag = sparse(n,1);
    sup_diag1 = sparse(n,1); % sup-diagonal +1
    sup_diag3 = sparse(n,1); % sup-diagonal +3

    % Diagonal elements
    main_diag(1) = 120 * x(1)^2 + 2 - 240 * x(1) * x(4) + 120 * x(4)^2; 
    main_diag(2) = 200 - 48 * x(2) * x(3) + 48 * x(3)^2 + 12 * x(2)^2; 
    main_diag(n-1) = 10 + 48 * x(n-2)^2 - 192 * x(n-2) * x(n-1) + 192 * x(n-1)^2; 
    main_diag(n) = 10 + 120 * x(n-3)^2 - 240 * x(n-3) * x(n) + 120 * x(n)^2;  
    
    % --- Special elements --- %
    sup_diag1(2) = 20; 
    sup_diag3(4) = -120 * x(1)^2 + 240 * x(1) * x(4) - 120 * x(4)^2;
    sup_diag1(3) = -24 * x(2)^2 + 96 * x(2) * x(3) - 96 * x(3)^2;  
    sup_diag1(n) = -10;

    % --- odd i = 3,5,...,n-3 --- %
    idx_odd = 3:2:n-3;
    main_diag(idx_odd) = 12 + 48 * x(idx_odd-1).^2 - 192 * x(idx_odd-1) .* x(idx_odd) + 312 * x(idx_odd).^2 - 240 * x(idx_odd) .* x(idx_odd+3) + 120 * x(idx_odd+3).^2;
    sup_diag1(idx_odd+1) = 10;
    sup_diag3(idx_odd+3) = -120 * x(idx_odd).^2 + 240 * x(idx_odd) .* x(idx_odd+3) - 120 * x(idx_odd+3).^2; 
    
    % --- even i = 4,6,...,n-2 --- %
    main_diag(idx_odd+1) = 210 + 120 * x(idx_odd-2).^2 - 240 * x(idx_odd-2) .* x(idx_odd+1) + 132 * x(idx_odd+1).^2 - 48 * x(idx_odd+1) .* x(idx_odd+2) + 48 * x(idx_odd+2).^2;
    sup_diag1(idx_odd+2) = -24 * x(idx_odd+1).^2 + 96 * x(idx_odd+1) .* x(idx_odd+2) - 96 * x(idx_odd+2).^2;
    sup_diag1(n-1) = -24 * x(n-2)^2 + 96 * x(n-2) * x(n-1) - 96 * x(n-1)^2; 
    
    sub_diag1 = [sup_diag1(2:n); 0];
    sub_diag3 = [sup_diag3(4:n); zeros(3,1)];
     
    % Final construction of the sparse hessian matrix
    H = spdiags([sub_diag3, sub_diag1, main_diag, sup_diag1, sup_diag3], ...
                [-3, -1, 0, 1, 3], n, n);
end
hess_powell = @(x,h) hess_powell_fun(x,h);

% Implementation of centered finite differences approximation for gradf and
% hessf:  
function g = gradf_finite_diff_powell_fun(x,h)
    n = length(x);
    g = sparse(n,1);
    g(1) = 40 * x(1)^3 + 2 * x(1) + 20 * x(2) - 120 * x(1)^2 * x(4) + 120 * x(1) * x(4)^2 - 40 * x(4)^3 + 40 * h(1)^2 * x(1);
    g(2) = 20 * x(1) + 200 * x(2) + 4 * x(2)^3 - 24 * x(2)^2 * x(3) + 48 * x(2) * x(3)^2 - 32 * x(3)^3 - 8 * h(2)^2 * x(3); 
    
    idx_odd = 3:2:n-3;
    
    % elements 3:2:n-3
    g(idx_odd) = 12 * x(idx_odd) + 10 * x(idx_odd+1) - 8 * x(idx_odd-1).^3 + 48 * x(idx_odd-1).^2 .* x(idx_odd) - ...
        96 * x(idx_odd-1) .* x(idx_odd).^2 + 104 * x(idx_odd).^3 - 120 * x(idx_odd).^2 .* x(idx_odd+3) + 120 * x(idx_odd) .* x(idx_odd+3).^2 - ...
        40 * x(idx_odd+3).^3 + 104 * h(idx_odd).^2 .* x(idx_odd) - 40 * h(idx_odd).^2 .* x(idx_odd+2); 

    % elements 4:2:n-2
    g(idx_odd+1) = 10 * x(idx_odd) + 210 * x(idx_odd+1) - 40 * x(idx_odd-2).^3 + ...
        120 * x(idx_odd-2).^2 .* x(idx_odd+1) - 120 * x(idx_odd-2) .* x(idx_odd+1).^2 + ...
        44 * x(idx_odd+1).^3 - 24 * x(idx_odd+1).^2 .* x(idx_odd+2) + 48 * x(idx_odd+1) .* x(idx_odd+2).^2 - 32 * x(idx_odd+2).^3 + ...
        - 40 * h(idx_odd+1).^2 .* x(idx_odd-2) + 40 * h(idx_odd+1).^2 .* x(idx_odd+1) - 16 * h(idx_odd+1).^2 .* x(idx_odd+2); 

    g(n-1) =  10 * x(n-1) - 10 * x(n) - 8 * x(n-2)^3 + 48 * x(n-2)^2 * x(n-1) - 96 * x(n-2) * x(n-1)^2 + 64 * x(n-1)^3 + 64 * h(n-1)^2 * x(n-1);
    g(n) = - 10 * x(n-1) + 10 * x(n) - 40 * x(n-3)^3 + 120 * x(n-3)^2 * x(n) - 120 * x(n-3) * x(n)^2 + 40 * x(n)^3 - 40 * h(n)^2 * x(n-3) + 40 * h(n)^2 * x(n);

end
gradf_finite_diff_powell = @(x,h) gradf_finite_diff_powell_fun(x,h);

function H = hess_finite_diff_powell_fun(x,h)
    n = length(x);
     
    % Pre-allocate vectors for diagonals
    main_diag = sparse(n,1);
    sup_diag1 = sparse(n,1); % sup-diagonal +1
    sup_diag3 = sparse(n,1); % sup-diagonal +3

    % Diagonal elements
    main_diag(1) = 120 * x(1)^2 + 2 - 240 * x(1) * x(4) + 120 * x(4)^2 + 20 * h(1)^2; 
    main_diag(2) = 200 - 48 * x(2) * x(3) + 48 * x(3)^2 + 12 * x(2)^2 + 2 * h(2)^2; 
    main_diag(n-1) = 10 + 48 * x(n-2)^2 - 192 * x(n-2) * x(n-1) + 192 * x(n-1)^2 + 2 * h(n-1)^2; 
    main_diag(n) = 10 + 120 * x(n-3)^2 - 240 * x(n-3) * x(n) + 120 * x(n)^2 + 20 * h(n)^2;  
    
    % --- Special elements --- %
    sup_diag1(2) = 20; 
    sup_diag3(4) = -120 * x(1)^2 + 240 * x(1) * x(4) - 120 * x(4)^2 - 20 * h(4)^2;
    sup_diag1(3) = -24 * x(2)^2 + 96 * x(2) * x(3) - 96 * x(3)^2 -16 * h(3)^2 - 48 * h(3) * x(3) + 24 * h(2) * x(2);  
    sup_diag1(n) = -10;

    % --- odd i = 3,5,...,n-3 --- %
    idx_odd = 3:2:n-3;
    main_diag(idx_odd) = 12 + 48 * x(idx_odd-1).^2 - 192 * x(idx_odd-1) .* x(idx_odd) + 312 * x(idx_odd).^2 - 240 * x(idx_odd) .* x(idx_odd+3) + 120 * x(idx_odd+3).^2 + 22 * h(idx_odd).^2;
    sup_diag1(idx_odd+1) = 10;
    sup_diag3(idx_odd+3) = -120 * x(idx_odd).^2 + 240 * x(idx_odd) .* x(idx_odd+3) - 120 * x(idx_odd+3).^2 - 20 * h(idx_odd+3).^2; 
    
    % --- even i = 4,6,...,n-2 --- %
    main_diag(idx_odd+1) = 210 + 120 * x(idx_odd-2).^2 - 240 * x(idx_odd-2) .* x(idx_odd+1) + 132 * x(idx_odd+1).^2 - 48 * x(idx_odd+1) .* x(idx_odd+2) + 48 * x(idx_odd+2).^2 + 22 * h(idx_odd+1).^2;
    sup_diag1(idx_odd+2) = -24 * x(idx_odd+1).^2 + 96 * x(idx_odd+1) .* x(idx_odd+2) - 96 * x(idx_odd+2).^2 - 48 * h(idx_odd+2).*x(idx_odd+2) + 24 * h(idx_odd+1) .* x(idx_odd+1)- 16 * h(idx_odd+2).^2;
    sup_diag1(n-1) = -24 * x(n-2)^2 + 96 * x(n-2) * x(n-1) - 96 * x(n-1)^2 - 16 * h(n-1)^2 - 48 * h(n-1) * x(n-1) + 24 * h(n-2) * x(n-2); 
    
    sub_diag1 = [sup_diag1(2:n); 0];
    sub_diag3 = [sup_diag3(4:n); zeros(3,1)];
     
    % Final construction of the sparse hessian matrix
    H = spdiags([sub_diag3, sub_diag1, main_diag, sup_diag1, sup_diag3], ...
                [-3, -1, 0, 1, 3], n, n);
end
hess_finite_diff_powell = @(x,h) hess_finite_diff_powell_fun(x,h);

function H = hess_finite_diff_powell_fun_h_var(x,h)
    n = length(x);
     
    % Pre-allocate vectors for diagonals
    main_diag = sparse(n,1);
    sup_diag1 = sparse(n,1); % sup-diagonal +1
    sup_diag3 = sparse(n,1); % sup-diagonal +3

    % Diagonal elements
    main_diag(1) = 120 * x(1)^2 + 2 - 240 * x(1) * x(4) + 120 * x(4)^2 + 20 * h(1)^2; 
    main_diag(2) = 200 - 48 * x(2) * x(3) + 48 * x(3)^2 + 12 * x(2)^2 + 2 * h(2)^2; 
    main_diag(n-1) = 10 + 48 * x(n-2)^2 - 192 * x(n-2) * x(n-1) + 192 * x(n-1)^2 + 2 * h(n-1)^2; 
    main_diag(n) = 10 + 120 * x(n-3)^2 - 240 * x(n-3) * x(n) + 120 * x(n)^2 + 20 * h(n)^2;  
    
    % --- Special elements --- %
    sup_diag1(2) = 20; 
    sup_diag3(4) = -120 * x(1)^2 + 240 * x(1) * x(4) - 120 * x(4)^2 - 40 * h(1)^2 - 40 * h(4)^2 + 60 * h(1) * h(4) + 120 * h(4) * (x(1) - x(4)) - 120 * h(1) * (x(1) - x(4)); 
    sup_diag1(3) = -24 * x(2)^2 + 96 * x(2) * x(3) - 96 * x(3)^2 -16 * h(3)^2 - 48 * h(3) * x(3) + 24 * h(2) * x(2);  
    sup_diag1(n) = -10;

    % --- odd i = 3,5,...,n-3 --- %
    idx_odd = 3:2:n-3;
    main_diag(idx_odd) = 12 + 48 * x(idx_odd-1).^2 - 192 * x(idx_odd-1) .* x(idx_odd) + 312 * x(idx_odd).^2 - 240 * x(idx_odd) .* x(idx_odd+3) + 120 * x(idx_odd+3).^2 + 22 * h(idx_odd).^2;
    sup_diag1(idx_odd+1) = 10;
    sup_diag3(idx_odd+3) = -120 * x(idx_odd).^2 + 240 * x(idx_odd) .* x(idx_odd+3) - 120 * x(idx_odd+3).^2 - 40 * h(idx_odd).^2 - 40 * h(idx_odd+3).^2 ...
        + 60 * h(idx_odd) .* h(idx_odd+3) + 120 * h(idx_odd+3) .* (x(idx_odd) - x(idx_odd+3)) - 120 * h(idx_odd) .* (x(idx_odd) - x(idx_odd+3)); 
    
    % --- even i = 4,6,...,n-2 --- %
    main_diag(idx_odd+1) = 210 + 120 * x(idx_odd-2).^2 - 240 * x(idx_odd-2) .* x(idx_odd+1) + 132 * x(idx_odd+1).^2 - 48 * x(idx_odd+1) .* x(idx_odd+2) + 48 * x(idx_odd+2).^2 + 22 * h(idx_odd+1).^2;
    sup_diag1(idx_odd+2) = -24 * x(idx_odd+1).^2 + 96 * x(idx_odd+1) .* x(idx_odd+2) - 96 * x(idx_odd+2).^2 ...
        -16 * h(idx_odd+2).^2 - 48 * h(idx_odd+2) .* x(idx_odd+2) + 24 * h(idx_odd+1) .* x(idx_odd+1);  
    
    sup_diag1(n-1) = -24 * x(n-2)^2 + 96 * x(n-2) * x(n-1) - 96 * x(n-1)^2 - 8 * h(n-2)^2 - 32 * h(n-1)^2 + 24 * h(n-2) * h(n-1) + 24 * x(n-2) * (2*h(n-1) - h(n-2)) + 48 * x(n-1) * (h(n-2) - 2 * h(n-1)); 
    
    sub_diag1 = [sup_diag1(2:n); 0];
    sub_diag3 = [sup_diag3(4:n); zeros(3,1)];
     
    % Final construction of the sparse hessian matrix
    H = spdiags([sub_diag3, sub_diag1, main_diag, sup_diag1, sup_diag3], ...
                [-3, -1, 0, 1, 3], n, n);
end
hess_finite_diff_powell_h_var = @(x,h) hess_finite_diff_powell_fun_h_var(x,h);


% function H = hess_finite_diff_powell_h_variable(grad, x, h)
%     n = length(x);
%     H = sparse(n,n);
%     gradx = grad(x);
% 
%     % select gradient in positions: 1, 6, 9, 14, 17...
%     jumps1 = repmat([5, 3], 1, ceil(n/2));
%     el1 = cumsum([1, jumps1]); % starts from position 1
%     el1(el1 > n) = []; 
%     x1 = x;
%     x1(el1) = x(el1) + h(el1);
% 
%     % select gradient in positions: 2, 5, 10, 13, 18...
%     jumps2 = repmat([3, 5], 1, ceil(n/2));
%     el2 = cumsum([2, jumps2]); % starts from position 2
%     el2(el2 > n) = []; 
%     x2 = x;
%     x2(el2) = x(el2) + h(el2);
% 
%     % select gradient in positions: 3, 8, 11, 16, 19...
%     el3 = cumsum([3, jumps1]); % starts from position 3
%     el3(el3 > n) = []; 
%     x3 = x;
%     x3(el3) = x(el3) + h(el3);
% 
%     % select gradient in positions:4, 7, 12, 15, 20...
%     el4 = cumsum([4, jumps2]); % starts from position 4
%     el4(el4 > n) = []; 
%     x4 = x;
%     x4(el4) = x(el4) + h(el4);
% 
%     col1 = [];
%     h1 =  h(el1);
%     h1_corr = [];
%     for i = 1:length(el1)
%         if i == 1
%             col1 = [col1, el1(1), el1(1)];
%             h1_corr = [h1_corr, h1(1), h1(1)];
%         elseif mod(i,2) == 0 
%             col1 = [col1, el1(i), el1(i-1), el1(i), el1(i), el1(i)];
%             h1_corr = [h1_corr, h1(i), h1(i-1), h1(i), h1(i), h1(i)];
%         else
%             col1 = [col1, el1(i), el1(i), el1(i)];
%             h1_corr = [h1_corr, h1(i), h1(i), h1(i)];
%         end
%     end
% 
%     col2 = [];
%     h2 =  h(el2);
%     h2_corr = [];
%     for i = 1:length(el2)
%         if i == 1
%             col2 = [col2, el2(i), el2(i), el2(i)];
%             h2_corr = [h2_corr, h2(i), h2(i), h2(i)];
%         elseif mod(i,2) == 0
%             col2 = [col2, el2(i), el2(i), el2(i)];
%             h2_corr = [h2_corr, h2(i), h2(i), h2(i)];
%         else
%             col2 = [col2, el2(i), el2(i-1), el2(i), el2(i), el2(i)];
%             h2_corr = [h2_corr, h2(i), h2(i-1), h2(i), h2(i), h2(i)];
%         end
%     end
% 
%     col3 = [];
%     h3 =  h(el3);
%     h3_corr = [];
%     for i = 1:length(el3)
%         if mod(i,2) == 0
%             col3 = [col3, el3(i), el3(i-1), el3(i), el3(i), el3(i)];
%             h3_corr = [h3_corr, h3(i), h3(i-1), h3(i), h3(i), h3(i)];
%         else
%             col3 = [col3, el3(i), el3(i), el3(i)];
%             h3_corr = [h3_corr, h3(i), h3(i), h3(i)];
%         end
%     end
% 
%     col4 = [];
%     h4 =  h(el4);
%     h4_corr = [];
%     for i = 1:length(el4)
%         if i == 1
%             col4 = [col4, el4(i)*ones(1,4)];
%             h4_corr = [h4_corr, h4(i)*ones(1,4)];
%         elseif mod(i,2) == 0
%             col4 = [col4, el4(i)*ones(1,3)];
%             h4_corr = [h4_corr, h4(i)*ones(1,3)];
%         else
%             col4 = [col4, el4(i), el4(i-1), el4(i)*ones(1,3)];
%             h4_corr = [h4_corr, h4(i), h4(i-1), h4(i)*ones(1,3)];
%         end
%     end
% 
%     grad1 = grad(x1);
%     grad2 = grad(x2);
%     grad3 = grad(x3);
%     grad4 = grad(x4);
% 
%     H(1,col1(1)) = (grad1(1) - gradx(1))./(h1_corr(1));
%     H(1,col2(1)) = (grad2(1) - gradx(1))./(h2_corr(1));
%     H(1,col4(1)) = (grad4(1) - gradx(1))./(h4_corr(1));
%     for i=2:n-2
%         H(i,col1(i)) = (grad1(i) - gradx(i))./(h1_corr(i));
%         H(i,col2(i)) = (grad2(i) - gradx(i))./(h2_corr(i));
%         H(i,col3(i)) = (grad3(i) - gradx(i))./(h3_corr(i));
%         H(i,col4(i)) = (grad4(i) - gradx(i))./(h4_corr(i));
%     end
%     H(end-1,col1(end)) = (grad1(end-1) - gradx(end-1))./(h1_corr(end));
%     H(end-1,col3(end)) = (grad3(end-1) - gradx(end-1))./(h3_corr(end-1));
%     H(end,col3(end)) = (grad3(end) - gradx(end))./(h3_corr(end));
%     H(end-1,col4(end)) = (grad4(end-1) - gradx(end-1))./(h4_corr(end));
% 
% 
%     H = 0.5 * (H + H');
% end 
% hess_finite_diff_powell_h_var = @(grad, x,h) hess_finite_diff_powell_h_variable(grad, x, h);

n = 1:dim; 
x3_powell = zeros(dim, 1); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(dim,1);

l_bound = x3_powell - 1;
u_bound = x3_powell + 1; 
M_ten_initial_points = l_bound + (u_bound - l_bound) .* rand(dim, 10);
starting_points = [x3_powell, M_ten_initial_points];

prinf_function(x3_opt, f3_powell)

rho = 0.9;
tolgrad = 10^(-3);
c1 = 10^(-4); 
btmax = 50;
kmax = 200;

flag_h = 0;
h = 0;

diff_pow = zeros(11,1);
time_pow = zeros(11, 1);
iter_pow = zeros(11, 1);
rate_pow = zeros(1,11); 
failure_pow = zeros(1,11);
for index=1:11
    tic;
    [xk, fk, gradfk_norm, k, ...
        xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f3_powell, grad_powell, hess_powell, kmax, ...
                                                            tolgrad, c1, btmax, dim, rho, h, flag_h);
    time_pow(index) = toc;
    diff_pow(index) = norm(xk - x3_opt,2);
    iter_pow(index) = k;
    rate_pow(index) = compute_exp_rate_conv_multi(xseqk);
    failure_pow(index) = failure;
    disp(gradfk_norm)
end
final_rate_exact_deriv_pow = mean(rate_pow, 'omitnan');

% Finite differences
btmax= 50; 
% parameter used for finite differences
h_vec = 10.^[-2, -4, -6, -8, -10, -12];

% Testing different values of h, fixed and varying from point to point 
flag_h2 = 1; % h must be adapted to each component of x

time1_finite_diff_pow = zeros(length(h_vec),11);
rate_h_pow = zeros(1, 11);
final_rate_pow = zeros(1, length(h_vec));
norm_err_conv_pow = zeros(11, length(h_vec));
iters_pow = zeros(11, length(h_vec));
failure1_pow = zeros(length(h_vec),11);

time2_finite_diff_pow = zeros(length(h_vec),11);
rate_h2_pow = zeros(1, 11);
final_rate2_pow = zeros(1,length(h_vec));
norm_err_conv2_pow = zeros(11, length(h_vec));
iters2_pow = zeros(11, length(h_vec));
failure2_pow = zeros(length(h_vec),11);

final_point_temporary_pow = zeros(dim,11);
final_point_pow = cell(1,11);
final_point_temporary2_pow = zeros(dim,11);
final_point2_pow = cell(1,11);

for i=1:length(h_vec)
    disp(i)
    for index=1:11
        disp(index)
        % tic;
        % [xk, fk, gradfk_norm, k, ...
        %     xseqk, btseq, failure] = modified_newton_method(starting_points(:,index), f3_powell, ...
        %                                                     gradf_finite_diff_powell, hess_finite_diff_powell, kmax,...
        %                                                     tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h);
        % hold on;
        % time1_finite_diff_pow(i, index) = toc;
        % iters_pow(index, i) = k;
        % norm_err_conv_pow(index,i) = norm(x3_opt - xk,2);
        % rate_h_pow(index) = compute_exp_rate_conv_multi(xseqk);
        % final_point_temporary_pow(:,index) = xk;
        % failure1_pow(index, i) = failure;
        % disp(rate_h_pow(index))

        tic;
        [xk2, fk2, gradfk_norm2, k2, ...
            xseqk2, btseq2, failure2] = modified_newton_method(starting_points(:,index), f3_powell, ...
                                                                gradf_finite_diff_powell, hess_finite_diff_powell_h_var, kmax,...
                                                                tolgrad, c1, btmax, dim, rho, h_vec(i), flag_h2);
        time2_finite_diff_pow(i, index) = toc;
        iters2_pow(index, i) = k2;
        norm_err_conv2_pow(index, i) = norm(x3_opt - xk2,2);
        rate_h2_pow(index) = compute_exp_rate_conv_multi(xseqk2);
        final_point_temporary2_pow(:,index) = xk2;
        failure2_pow(index, i) = failure2;
    end
    final_rate_pow(i) = mean(rate_h_pow, 'omitnan');
    final_rate2_pow(i) = mean(rate_h2_pow, 'omitnan');
    final_point_pow{i} = final_point_temporary_pow;
    final_point2_pow{i} = final_point_temporary2_pow;
end

filename = 'output.mat';
%base_vars = {'time_pow', 'diff_pow', 'iter_pow', 'final_rate_exact_deriv_pow', 'failure_pow'};
%base_vars = {'time1_finite_diff_pow', 'norm_err_conv_pow', 'iters_pow', 'final_rate_pow', 'failure1_pow'};
base_vars = {'time2_finite_diff_pow', 'norm_err_conv2_pow', 'iters2_pow', 'final_rate2_pow', 'failure2_pow'};

for i = 1:length(base_vars)
    base_name = base_vars{i};
    new_name = sprintf('%s_%d', base_name, dim);  % es. time_pow_1000

    val = eval(base_name);

    eval([new_name ' = val;']);

    if exist(filename, 'file')
        save(filename, new_name, '-append');
    else
        save(filename, new_name);
    end
end