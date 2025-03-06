% Function to visualize Nelder Mead optimization on a 2D Rosenbrock function.
% It shows the trajectories of simplex's through iterations and initial points
% A background color map that represents the "depth" of the objective
% function helps to identify the region where the method should converge.

function nelderMead_picture2D(f, x_interval, y_interval, x, initial_point, type_of_move)
[X, Y] = meshgrid(x_interval, y_interval);   
Z = f(X, Y);
figure;
imagesc(x_interval, y_interval, Z);
set(gca, 'YDir', 'normal'); 
colorbar;
set(gca, 'ColorScale', 'log');  
hold on;
contour(X, Y, Z, 50, 'LineColor', 'k');  
for i = 1:size(x, 1)
    scatter(x(i, 1), x(i, 2), 35, type_of_move(i, :), 'filled');
end
plot(x(:, 1), x(:, 2), 'k-', 'LineWidth', 1.2);
plot(1, 1, 'kp', 'MarkerFaceColor', 'k', 'MarkerSize', 13.5, 'MarkerFaceColor', 'w'); 
plot(initial_point(1), initial_point(2), 'kp', 'MarkerSize', 13.5, 'MarkerFaceColor', 'r'); 
hold off;
xlabel('x');
ylabel('y');
title('Rosenbrock function and iteration with Nelder Mead');
end

