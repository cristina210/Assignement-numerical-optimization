% Tuning of parameters rho, sigma, gamma, chi used in Nelder Mead method.
% 3 different values of each parameters are compared in each dimensione =
% 10,25,50. In each section a different function of exercise 3 is take into
% account together with the initial starting point suggested from the pdf
% file.
% The comparison between parameters is accomplished through the study of
% two quantities: number of iteration to converge and error of convergence
% (distance from optimal point)
% The choice for the parameter is done take into consideration these
% two quantities: best parameter shoud minimize both of them.
% It is exploited function "comparePar2()"

%% Tuning in Chained Rosenbrock
clc
clear all
close all

% Value of parameters compared
l = 5; 
rho_vec = [0.25, 0.5, 1, 1.35, 1.75];
sigma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
gamma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
chi_vec = [1.1, 1.5, 2, 2.5, 3];

% Comparison in different dimension
vec_dim = [10,25,50];

for j = 1:length(vec_dim)
    
    disp("Dimension:")
    dim = vec_dim(j);
    disp(dim)
    
    % Rosenbrock function
    % Information of function and initial point (suggested from PDF)
    f= @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

    initial_point = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

    x_opt = ones(1,dim);

    % Nelder Mead
    [pos_bestRho, pos_bestSigma, pos_bestGamma, pos_bestChi] = comparePar2(dim, f, initial_point,x_opt, l, rho_vec,sigma_vec, gamma_vec,chi_vec);
    disp("Best configuration:")
    disp("value of rho")
    disp(rho_vec(pos_bestRho))
    disp("value of sigma")
    disp(sigma_vec(pos_bestSigma))
    disp("value of gamma")
    disp(gamma_vec(pos_bestGamma))
    disp("value of chi")
    disp(chi_vec(pos_bestChi))
end  

%% Tuning in Wood's function
clc
clear all
close all

% Value of parameters compared
l = 5; 
rho_vec = [0.25, 0.5, 1, 1.35, 1.75];
sigma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
gamma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
chi_vec = [1.1, 1.5, 2, 2.5, 3];

% Comparison in different dimension
vec_dim = [10,25,50];

for j = 1:length(vec_dim)
    dim = vec_dim(j);
    disp(dim)
    
    % Wood function
    % Information of function and initial point (suggested from PDF)
    f= @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

    n = 1:dim; 
    initial_point = zeros(1,dim); 
    initial_point(mod(n,2) == 1 & n <= 4) = -3;
    initial_point(mod(n,2) == 1 & n > 4) = -2; 
    initial_point(mod(n,2) == 0 & n <= 4) = -1;
    initial_point(mod(n,2) == 0 & n > 4) = 0; 

    x_opt = ones(1,dim);

    % Nelder Mead
    [pos_bestRho, pos_bestSigma, pos_bestGamma, pos_bestChi] = comparePar2(dim, f, initial_point, x_opt, l, rho_vec,sigma_vec, gamma_vec,chi_vec);
    disp("Best configuration:")
    disp("value of rho")
    disp(rho_vec(pos_bestRho))
    disp("value of sigma")
    disp(sigma_vec(pos_bestSigma))
    disp("value of gamma")
    disp(gamma_vec(pos_bestGamma))
    disp("value of chi")
    disp(chi_vec(pos_bestChi))
end 


%% Tuning in Powell function
clc
clear all
close all

% Value of parameters compared
l = 5; 
rho_vec = [0.25, 0.5, 1, 1.35, 1.75];
sigma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
gamma_vec = [0.1, 0.25, 0.5, 0.75, 0.9];
chi_vec = [1.1, 1.5, 2, 2.5, 3];

% Comparison in different dimension
vec_dim = [10,25,50];


for j = 1:length(vec_dim)
    
    disp("Dimension:")
    dim = vec_dim(j);
    disp(dim)
    
    % Powell function
    % Information of function and initial point (suggested from PDF)
    f = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

    n = 1:dim; 
    initial_point = zeros(1, dim); 
    initial_point(mod(n,4) == 1) = 3;   
    initial_point(mod(n,4) == 2) = -1; 
    initial_point(mod(n,4) == 3) = 0;  
    initial_point(mod(n,4) == 0) = 1; 

    x_opt = zeros(1,dim);

    % Nelder Mead
    [pos_bestRho, pos_bestSigma, pos_bestGamma, pos_bestChi] = comparePar2(dim, f, initial_point,x_opt, l, rho_vec,sigma_vec, gamma_vec,chi_vec);
    disp("Best configuration:")
    disp("value of rho")
    disp(rho_vec(pos_bestRho))
    disp("value of sigma")
    disp(sigma_vec(pos_bestSigma))
    disp("value of gamma")
    disp(gamma_vec(pos_bestGamma))
    disp("value of chi")
    disp(chi_vec(pos_bestChi))
end  

