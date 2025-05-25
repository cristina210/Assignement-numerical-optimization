function rate = compute_exp_rate_conv_multi(xseq)

% function rate = compute_exp_rate_conv_multi(xseq)
%
% Function that estimates the experimental rate of convergence (ERC) 
% of an iterative optimization method, based on the last three iterates.
%
% INPUT:
% xseq = n-by-m matrix where each column represents the iterate x_k of the 
%        optimization method. 
%
% OUTPUT:
% rate = estimate of the experimental rate of convergence;

num = xseq(:,end) - xseq(:,end-1);
denom = xseq(:,end-1) - xseq(:,end-2);
rate = log(sum(num.^2))/log(sum(denom.^2));

end
