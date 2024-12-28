function plot_areaerrorbar(data, x_axis)

% Default options
options.handle     = figure(1);
%         options.color_area = [128 193 219]./255;    % Blue theme
options.color_line = [ 52 148 186]./255;
options.color_area = [243 169 114]./255;    % Orange theme
%options.color_line = [236 112  22]./255;
options.alpha      = 0.5;
options.line_width = 2;

options.x_axis = x_axis';

% Computing the mean and standard deviation of the data matrix
data_mean = mean(data,1);
error  = max(data)-min(data);

% Plotting the result
figure(options.handle);
x_vector = [options.x_axis', fliplr(options.x_axis')];
patch = fill(x_vector, [data_mean+error/2,fliplr(data_mean-error/2)], options.color_area);
set(patch, 'edgecolor', 'none');
set(patch, 'FaceAlpha', options.alpha);
hold on;
plot(options.x_axis, data_mean, 'color', options.color_line, ...
    'LineWidth', options.line_width);
hold off;

end