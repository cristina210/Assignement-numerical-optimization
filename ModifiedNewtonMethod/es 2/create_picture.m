function create_picture(matrix, rho_vec, c1_vec, point, val_y, value_studied)
%
% function create_picture(matrix, rho_vec, c1_vec, point)
% Function that create a bar plot showing how the number of iterations
% needed to converge vary base on the different combination of rho and c1
% values
% 
% INPUTS:
% matrix = matrix with length(rho_vec) rows and length(c1_vec) columns
% whose entries are the iterations required for each values combination
% rho_vec = vector containing all the rho values
% c1_vec = vector containing all the c1 values
% point = term which specifies which is the starting point (1,2)

    figure;

    data = matrix';
    b = bar(data, 'grouped');  
    colors = lines(length(rho_vec));
    for i = 1:length(rho_vec)
        b(i).FaceColor = colors(i, :);
    end

    xticks(1:length(c1_vec));
    xticklabels(arrayfun(@(x) sprintf('%.0e', x), c1_vec, 'UniformOutput', false));
    xlabel(val_y, 'FontSize', 14);
    ylabel(value_studied, 'FontSize', 14);
    legend(arrayfun(@(x) sprintf('\\rho = %.1f', x), rho_vec, 'UniformOutput', false), ...
           'Location', 'BestOutside');
    title(sprintf('%s vs %s from point %d', value_studied, val_y, point), 'FontSize', 16);
    grid on;
    
end
