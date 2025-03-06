% Point 3, analysis of Nelder-Mead method 
% Script for Evaluating the Nelder-Mead method.
% This script implements and tests the Nelder-Mead algorithm
% on three multidimensional benchmark functions:
%  - Chained Rosenbrock's Function
%  - Wood's Function
%  - Powell's Function
%
% For each function, simulations are performed to analyze the convergence, 
% computational cost, stagnation and empirical rates of convergence of the method in
% various dimensions (10,25,50) and multiple initial points (one suggested from 
% PDF and 10 points from the ipercube).
% It exploits NelderMead_for10Points() function.

%% Chained Rosenbrook per dim = 10 %%
clc
clear all
close all
rng(min(343341,343428))  
dim = 10;
disp("dimensione:")
disp(dim)
rho = 1;
sigma = 0.1;
chi = 1.5;
gamma = 0.75;

f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f1_ros,x1_rosenbrock,x1_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from 10 points from the ipercube:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)  

vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % contains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x1_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x1_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Chained Rosenbrook per dim = 25 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 25;
disp("dimensione:")
disp(dim)
rho = 1;
sigma = 0.1;
chi = 2;
gamma = 0.5;

f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f1_ros,x1_rosenbrock,x1_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x1_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x1_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Chained Rosenbrook per dim = 50 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 50;
disp("dimensione:")
disp(dim)
rho = 0.5;
sigma = 0.1;
chi = 2.5;
gamma = 0.9;

f1_ros = @(x) sum(arrayfun(@(i) 100*(x(i)^2 - x(i+1))^2 + (x(i) - 1)^2, 1:length(x)-1));

x1_rosenbrock = arrayfun(@(i) -1.2*(mod(i,2)==1) + 1.0*(mod(i,2)==0), 1:dim);

x1_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f1_ros,x1_rosenbrock,x1_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x1_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x1_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Wood function per dim = 10 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 10;
disp("dimensione:")
disp(dim)
rho = 1;
sigma = 0.5;
chi = 2.5;
gamma = 0.75;

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;  

x2_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points, x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f2_wood,x2_wood,x2_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x2_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x2_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end 
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Wood function per dim = 25 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 25;
disp("dimensione:")
disp(dim)
rho = 1;
sigma = 0.1;
chi = 2;
gamma = 0.9;

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;   
x2_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f2_wood,x2_wood,x2_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x2_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x2_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end 
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end


%% Wood function per dim = 50 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 50;
disp("dimensione:")
disp(dim) 
rho = 1;
sigma = 0.1;
chi = 2;
gamma = 0.9;

f2_wood = @(x) sum(arrayfun(@(j) ...
    100*(x(2*j-1)^2 - x(2*j))^2 + (x(2*j-1) - 1)^2 + ...
    90*(x(2*j+1)^2 - x(2*j+2))^2 + (x(2*j+1) - 1)^2 + ...
    10*(x(2*j) + x(2*j+2) - 2)^2 + (x(2*j) - x(2*j+2))^2 / 10, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x2_wood = zeros(1,dim); 
x2_wood(mod(n,2) == 1 & n <= 4) = -3;
x2_wood(mod(n,2) == 1 & n > 4) = -2; 
x2_wood(mod(n,2) == 0 & n <= 4) = -1;
x2_wood(mod(n,2) == 0 & n > 4) = 0;   
x2_opt = ones(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points,lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f2_wood,x2_wood,x2_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)


% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x2_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x2_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Powell function per dim = 10 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 10;
disp("dimensione:")
disp(dim)
rho = 1.3;
sigma = 0.1;
chi = 2;
gamma = 0.75;

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f3_powell,x3_powell,x3_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x3_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x3_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Powell function per dim = 25 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 25;
disp("dimensione:")
disp(dim)
rho = 1.35;
sigma = 0.1;
chi = 1.1;
gamma = 0.75;

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f3_powell,x3_powell,x3_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x3_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x3_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end

%% Powell function per dim = 50 %%
clc
clear all
close all
rng(min(343341,343428))
dim = 50;
disp("dimensione:")
disp(dim)
rho = 1;
sigma = 0.1;
chi = 1.1;
gamma = 0.9;

f3_powell = @(x) sum(arrayfun(@(j) ...
    (x(2*j-1) + 10*x(2*j))^2 + 5*(x(2*j+1) - x(2*j+2))^2 + ...
    (x(2*j) - 2*x(2*j+1))^4 + 10*(x(2*j-1) - x(2*j+2))^4, ...
    1:(length(x)-2)/2));

n = 1:dim; 
x3_powell = zeros(1, dim); 
x3_powell(mod(n,4) == 1) = 3;   
x3_powell(mod(n,4) == 2) = -1; 
x3_powell(mod(n,4) == 3) = 0;  
x3_powell(mod(n,4) == 0) = 1;  
x3_opt = zeros(1,dim);

[vec_time, k1, x_1, k_10_points,x_10_points, lista_rates, lista_err] = NelderMead_for_10Points_2(dim,f3_powell,x3_powell,x3_opt,rho,sigma,chi,gamma);

disp("Convergence point from initial point suggested from pdf file:")
disp(x_1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k1)
disp("Number iteration before convergence from initial point suggested from pdf file:")
disp(k_10_points)
disp("Convergence points from 10 initial points from the ipercube and flag:")
disp(x_10_points)
disp("Computational cost for each of the 1 + 10 points")
disp(vec_time)

% for latex output
vec_dist = zeros(1,11); % contains the distance from optimum of the convergence point with the 11 fixed initial point
type_conv = zeros(1,11);   % conitains the type of convergence (no convergence, convergence with 1^/2^ stopping criteria)
vec_dist(1,1) = norm(x_1(1:dim) - x3_opt);
type_conv(1,1) = x_1(dim + 1);
for i=1:11
    if i ~= 1
        dist_i = norm(x_10_points(i-1,1:dim) - x3_opt);
        vec_dist(1,i) = dist_i;
        type_conv(1,i) = x_10_points(i-1, dim + 1);
    end
end
disp("Distance from optimum for each of the 11 points")
disp(vec_dist)
for i = 1:length(lista_rates)
    fprintf("Last 5 experimental rates of convergence of point: %d\n", i);
    vec = lista_rates{i};  
    disp(vec(max(1, end-4):end));
end