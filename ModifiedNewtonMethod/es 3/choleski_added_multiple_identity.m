function L = choleski_added_multiple_identity(hessfk, beta, n)

% function L = choleski_added_multiple_identity(hessfk, beta, n)
% 
% Function that modifies a possibly indefinite Hessian matrix by adding a 
% scaled identity matrix to ensure positive definiteness, and computes its 
% Cholesky factorization.
% 
% INPUTS:
% hessfk = symmetric Hessian matrix (possibly indefinite) of size n-by-n;
% beta = small positive scalar used to ensure the perturbed matrix becomes
%        positive definite;
% n = dimension of the matrix hessfk;
% 
% OUTPUT:
% L = lower triangular matrix such that L * L' = Bk, where Bk is a 
%     positive definite modification of hessfk.
%
% DESCRIPTION:
% The function iteratively adds tau_k * I to hessfk, where tau_k is
% initialized based on the smallest diagonal entry of hessfk. If hessfk is 
% already positive definite, tau_k starts at zero. Otherwise, a positive 
% shift is applied.
% At each iteration, a Cholesky decomposition is attempted. If it fails
% (i.e., the matrix is not positive definite), tau_k is increased (up to a 
% maximum of 1e5) and the process is repeated, up to kmax iterations.
% If a positive definite matrix is found, the function returns its 
% Cholesky factor. Otherwise, it raises an error.

kmax = 500;

% Find minimun element of the diagonal
min_diag_hess = min(diag(hessfk));

% Initialization of tauk
if min_diag_hess > 0
    tauk = 0;
else
    tauk = -min_diag_hess + beta;
end

flag_pos_def = 1;
k = 0;
I = speye(n);

while k < kmax && flag_pos_def ~= 0
    % Bk is created as a sparse matrix
    Bk = hessfk + tauk * I;

    % Cholesky factoritazion to check if the matrix is positive definite
    [L, flag_pos_def] = chol(Bk, 'lower'); 

    if flag_pos_def == 0
        break;
    else
        % tauk = max(10 * tauk, beta);
        % tauk = max(2 * tauk, beta);  
        tauk = min(max(2 * tauk, beta), 1e5);  
    end

    k = k + 1;
end

if tauk > 1e6
    error('Tauk too big.')
elseif k == kmax && flag_pos_def == 1
    error('Bk NOT positive definite')
end

end