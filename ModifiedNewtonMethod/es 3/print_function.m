function print_function(x_opt, f)
    i = 1;  % index 1^st variable to vary
    j = 2;  % index 2^nd variable to vary

    % Meshgrid in the neighbourhood of x_opt 
    delta = 1;
    step = 0.1;
    x1 = x_opt(i);
    x2 = x_opt(j);
    [xi, xj] = meshgrid(x1 - delta : step : x1 + delta, ...
                        x2 - delta : step : x2 + delta);

    Z = zeros(size(xi));

    % Evaluation of f varying xi, xj
    for row = 1:size(xi, 1)
        for col = 1:size(xi, 2)
            x = x_opt;
            x(i) = xi(row, col);
            x(j) = xj(row, col);
            Z(row, col) = f(x);
        end
    end

    % Plot of the neighbourhood of x_opt
    figure;
    surf(xi, xj, Z);
    xlabel(sprintf('x_%d', i));
    ylabel(sprintf('x_%d', j));
    zlabel('f(x)');
    title(sprintf('Surrounding of x_%d, x_%d', i, j));
    shading interp;
    view(135,30);
    hold on;

    f_min = f(x_opt);
    plot3(x_opt(1), x_opt(2), f_min, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    legend('f(x)', 'theoric minimum point');
    hold off;

end
