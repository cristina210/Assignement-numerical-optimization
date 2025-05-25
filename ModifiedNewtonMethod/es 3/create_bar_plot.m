function create_bar_plot(vec1, vec2, vec3, flag, flag_fun)

% function create_bar_plot(vec1, vec2, vec3, flag)
% 
% This function creates a grouped bar plot to compare three vectors of 
% performance metrics (e.g., time, error, iterations, or convergence rate) 
% across different problem dimensions (e.g., 1e3, 1e4, 1e5). 
% If the vectors contain only one element each, the x-axis will represent 
% the dimensions directly. Otherwise, the x-axis will show values of the 
% step size (h), and the bars will be grouped by step size, with colors 
% indicating different dimensions.
%
% INPUTS:
%   vec1 - Vector of metrics for dimension 1e3
%   vec2 - Vector of metrics for dimension 1e4
%   vec3 - Vector of metrics for dimension 1e5
%   flag - A string indicating the type of metric being plotted 
%          ('time', 'diff', 'iter', or other)
    
    n = length(vec1);
    rates = [vec1(:), vec2(:), vec3(:)];  % Each column is a dimension

    figure;

    if n == 1
        for i = 1:3
            bar(i, rates(1, i), 'FaceColor', get_color(i));
            hold on;
        end
        xlim([0.5, 3.5]);
        xticks(1:3);
        xticklabels({'dim = 1e3', 'dim = 1e4', 'dim = 1e5'});
        xlabel('Dimension');
    else
        hb = bar(rates); 
        colors = [0.8500 0.3250 0.0980;   
                  0.4660 0.6740 0.1880;   
                  0.0000 0.4470 0.7410];  
        for i = 1:3
            hb(i).FaceColor = colors(i,:);
        end

        x_labels = arrayfun(@(k) sprintf('10^{-%d}', 2*k), 1:n, 'UniformOutput', false);
        xticks(1:n);
        xticklabels(x_labels);
        xlabel('h');
    end

    ylim([0, max(rates(:)) + 0.01]);

    switch flag
        case 'time' 
            ylabel('Average time');
            metric_title = 'Comparison of average time';
        case 'diff'
            ylabel('Average difference');
            metric_title = 'Comparison of average difference';
        case 'iter'
            ylabel('Average iterations');
            metric_title = 'Comparison of average iterations';
        otherwise
            ylabel('Average rate of convergence');
            metric_title = 'Comparison of average rate of convergence';
    end

    problem_name = '';
    switch flag_fun
        case 'ros'
            problem_name = ' (Rosenbrock Chained)';
        case 'wood'
            problem_name = ' (Wood Chained)';
        case 'pow'
            problem_name = ' (Powell Chained)';
    end

    title([metric_title, problem_name]);


    legend({'dim = 1e3', 'dim = 1e4', 'dim = 1e5'}, 'Location', 'best');
    grid on;
end

function c = get_color(i)

% function get_color(i)
% Helper function to return a predefined color for each dimension index.
%
% INPUTS:
%   i - Index of the vector (1, 2, or 3)
%
% OUTPUTS:
%   c - RGB color triplet
    color_map = [0.8500 0.3250 0.0980;  
                 0.4660 0.6740 0.1880;   
                 0.0000 0.4470 0.7410];  
    c = color_map(i, :);
end
