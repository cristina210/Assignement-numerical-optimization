function L = choleski_added_multiple_identity(hessfk, beta, n)
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