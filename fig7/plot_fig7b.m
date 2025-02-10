% Clear stuff
clear all; close all;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Load big PTC scan data
load('../data_files/fig7_data.mat');

%--------------------------------%
%     Coordinates for 'Hole'     %
%--------------------------------%
% G-direction
intersection.theta_old = 0.3176;
intersection.A_perturb = 0.5576;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
% Green     (#2ca02c) = [ 44, 160,  44] ./ 255
% Chartreus (#bcbd22) = [188, 189,  34] ./ 255
% Yellow    (#fafa2a) = [250, 250,  42] ./ 255
% Orange    (#ff7f0e) = [255, 126,  14] ./ 255
% Red       (#d62728) = [214,  39,  40] ./ 255
% Pink      (#e38ab7) = [227, 138, 183] ./ 255
% Purple    (#9467bd) = [148, 103, 189] ./ 255
% Cyan      (#1bc3cc) = [ 27, 195, 204 ./ 255

% Plot colours
plot_colours = {[188, 189,  34] ./ 255;
                [250, 250,  42] ./ 255;
                [255, 126,  14] ./ 255;
                [214,  39,  40] ./ 255;
                [227, 138, 183] ./ 255;
                [148, 103, 189] ./ 255;
                [ 27, 195, 204] ./ 255};

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Colourmap     %
%-------------------%
% Plotting colours
colours = colororder();

n_colours_in = 2048;
n_colours_out = 1024;

% Get colour map
colour_map = parula(n_colours_in);
% Create a linear ramp the size of the colormap we actually want
t = linspace(0,1,n_colours_out);
% Apply whatever transform you like to the ramp
t2 = t .^ 2;
% Use that to scale the big linear colormap into the small stretched one.
colour_map_transformed = colour_map(1+floor((n_colours_in-1)*t2'),:); 

colormap(colour_map_transformed);

%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(7); clf;
fig.Name = 'PTC Scans';
ax = gca();

% Axis dimensions
width = 8.0;
height = 6.4;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%-------------------------------------------------%
%     Plot: Stable Manifold Intersection Pole     %
%-------------------------------------------------%
plot3(ax, [intersection.theta_old, intersection.theta_old], ...
     [intersection.A_perturb, intersection.A_perturb], ...
     [-5, 5], ...
     Color='k', LineWidth=2.5, LineStyle='-');

%--------------------------%
%     Plot: PTC Curves     %
%--------------------------%
% Linewidth
lw = 1.5;

plot_A_perturb = [0.05, 0.1, 0.15, 0.55, 1.0, 1.5, 2.0];

% Find plotting indices
plot_idx = zeros(length(plot_A_perturb), 1);
for i = 1 : length(plot_A_perturb)
  plot_idx(i) = find(round(A_perturb, 3) == plot_A_perturb(i));
end

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % A_perturb plot data
  A_plot_lt1 = A_perturb(idx) * ones(length(theta_old_lt1{idx}));
  A_plot_gt1 = A_perturb(idx) * ones(length(theta_old_gt1{idx}));

  % Plot
  plot3(ax, theta_old_lt1{idx}, A_plot_lt1, theta_new_lt1{idx}, Color=plot_colours{i}, LineStyle='-');
  plot3(ax, theta_old_gt1{idx}-1, A_plot_gt1, theta_new_gt1{idx}, Color=plot_colours{i}, LineStyle='-');
end

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = pad_data(data_hole_lt1, 0, 'lt1');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = pad_data(data_hole_gt1, 0, 'gt1');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: Before hole
[X, Y, Z] = pad_data(data_before_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: After hole
[X, Y, Z] = pad_data(data_after_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

%---------------------------%
%     Surface: Settings     %
%---------------------------%
% Shading of surface
shading(ax, 'interp');

% Colorbar
% cbar = colorbar();

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 2.0];
ax.ZAxis.Limits = [0.0, 2.5];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 1.0 : 2.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 0.5 : 2.0;

% Z-Axis
ax.ZAxis.TickValues = 0.0 : 0.5 : 3.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 0.0 : 0.25 : 3.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% % Axis labels
% xlabel(ax, '$\theta_{\mathrm{o}}$');
% ylabel(ax, '$A_{\mathrm{p}}$');
% zlabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all axis labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};
ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
view(135, 15);
filename_out = '../fig7b_G_PTC_surface_2.png';

exportgraphics(fig, filename_out, ContentType='image', Resolution=750);

%-------------------------------------------------------------------------%
%%                           Data Functions                              %%
%-------------------------------------------------------------------------%
function [theta_old, A_perturb, theta_new] = pad_data(data_in, theta_new_modifier, theta_old_lt1_gt1)
  % [theta_old, theta_new, A_perturb] = pad_data(theta_old_in, theta_new_in, A_perturb_in)
  %
  % Read the data and pad it to the max length data array
  
  % Read data
  theta_old_read = data_in.theta_old;
  theta_new_read = data_in.theta_new;
  A_perturb_read = data_in.A_perturb;

  % Output data
  theta_old = [];
  theta_new = [];
  A_perturb = [];

  % Get max length
  lengths = [];
  max_theta_old = [];
  min_theta_old = [];
  for i = 1 : length(theta_old_read)
    % Read temp data
    theta_old_temp = theta_old_read{i};

    max_theta_old = [max_theta_old, max(theta_old_temp)];
    min_theta_old = [min_theta_old, min(theta_old_temp)];

    % Get max length
    lengths = [lengths, length(theta_old_temp)];
  end
  [~, max_idx] = max(lengths);

  % Get max length theta_old
  theta_old_max = theta_old_read{max_idx};

  if strcmp(theta_old_lt1_gt1, 'lt1')
    % Find theta_old_read with max theta_old value, interpolate so is the
    % same length as theta_old_max
    [~, max_idx] = min(min_theta_old);

    % Interpolate date
    theta_old_max = theta_old_read{max_idx};
  end

  if strcmp(theta_old_lt1_gt1, 'gt1')
    % Find theta_old_read with max theta_old value, interpolate so is the
    % same length as theta_old_max
    [~, max_idx] = max(max_theta_old);

    % Interpolate date
    theta_old_max = theta_old_read{max_idx};
  end

  % Cycle through all data arrays and interpolate theta_new
  for i = 1 : length(theta_new_read)

    % Read temp data
    theta_old_temp = theta_old_read{i};
    theta_new_temp = theta_new_read{i};
    A_perturb_temp = A_perturb_read(i);

    % Get unique indices
    [~, unique_idx] = unique(theta_old_temp);
    theta_old_temp = theta_old_temp(unique_idx);
    theta_new_temp = theta_new_temp(unique_idx);

    % Check if i != max_idx
    if i ~= max_idx
      % Interpolate data
      theta_new_interp = interp1(theta_old_temp, theta_new_temp, theta_old_max);
    else
      % Just append that array
      theta_new_interp = theta_new_read{i};
    end
    
    theta_old_interp = theta_old_max;

    % If doing data theta_old < 1, set all theta_old values < the min value
    % of that data set to the min value
    if strcmp(theta_old_lt1_gt1, 'lt1')
      % Find min value
      [min_val, min_idx] = min(theta_old_read{i});
      % Find all values less than min val
      lt1_idx = theta_old_interp < min_val;

      % % Find max theta_old value
      % [max_val, max_idx] = max(theta_old_read{i});
      % % Find all values greater than max val
      % gt1_idx = theta_old_interp > max_val;

      % Set to the value
      theta_old_interp(lt1_idx) = min_val;

      % Find NaNs
      is_theta_new_nan = isnan(theta_new_interp);

      % Set NaNs of lt1_idx to min val
      theta_new_interp(is_theta_new_nan(lt1_idx)) = theta_new_read{i}(min_idx);

      % Set NaNs of gt1_idx to max val
      % theta_new_interp(is_theta_new_nan(gt1_idx)) = theta_new_read{i}(max_idx);
    end

    if strcmp(theta_old_lt1_gt1, 'gt1')
      % Find min value
      [max_val, max_idx] = max(theta_old_read{i});

      % Find all values less than this
      gt1_idx = theta_old_interp > max_val;

      % Set to the value
      theta_old_interp(gt1_idx) = max_val;

      % Set all NaN values to the first theta_new value
      theta_new_NaN = theta_new_read{i}(max_idx);
      % Find NaNs
      theta_new_interp(isnan(theta_new_interp)) = theta_new_NaN;
    end
    
    % % Append the first value
    % theta_new_interp = [theta_new_1, theta_new_interp];
    % theta_old_interp = [theta_old_1, theta_old_max];

    % Append other stuff
    theta_new = [theta_new; theta_new_interp];
    theta_old = [theta_old; theta_old_interp];
    A_perturb = [A_perturb; A_perturb_read(i) * ones(1, length(theta_new_interp))];
  end

  % Add modifier
  theta_new = theta_new + theta_new_modifier;

end
