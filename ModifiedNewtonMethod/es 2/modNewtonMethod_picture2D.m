function modNewtonMethod_picture2D(f, x_interval, y_interval, xseq1, xseq2, initial_point1, initial_point2)
% 
% function modNewtonMethod_picture2D(f, x_interval, y_interval, xseq1, xseq2, initial_point1, initial_point2)
% 
% Function that take the picture whit color and countout map of the
% Rosenbrock function, showing the steps needed to the modified Newton
% method to converge to the minimum point startinf from the two different
% initial point given.
%
% INPUTS:
% f = function handle that describes a function R^n->R;
% x_interval = x coords in which the function is evaluated;
% y_interval = y coords in which the function is evaluated;
% xseq1 = sequences of points find each step until the method converges
% from initial_point1;
% xseq2 = sequences of points find each step until the method converges
% from initial_point2;
% initial_point1 = 1^st starting point
% initial_point2 = 2^nd starting point
% 

[X, Y] = meshgrid(x_interval, y_interval);  
Z = f(X, Y); 
figure;
imagesc(x_interval, y_interval, Z);
set(gca, 'YDir', 'normal'); 
colorbar;
set(gca, 'ColorScale', 'log');  
hold on;
contour(X, Y, Z, 50, 'LineColor', 'k');
plot(initial_point1(1), initial_point1(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r'); 
plot(initial_point2(1), initial_point2(2), 'cp', 'MarkerSize', 12, 'MarkerFaceColor', 'c'); 
plot(xseq1(1, :), xseq1(2, :), 'ro-', 'LineWidth', 1.3, 'MarkerSize', 3.5)
plot(xseq2(1, :), xseq2(2, :), 'co-', 'LineWidth', 1.3, 'MarkerSize', 3.5)
plot(1, 1,  'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 5); 
hold off;
xlabel('x');
ylabel('y');
title('Rosenbrock function - Color and Contour Map');

end

