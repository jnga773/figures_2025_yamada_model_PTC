% Clear stuff
clear all;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
mat_file = '../data_files/fig9_data.mat';
% Load data from .mat
load(mat_file);

%--------------------------------%
%     Coordinates for 'Hole'     %
%--------------------------------%
% I-Direction
intersection.theta_old = 0.503722;
intersection.A_perturb = 4.142870;

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
t2 = t .^ 0.9;
% Use that to scale the big linear colormap into the small stretched one.
colour_map_transformed = colour_map(1+floor((n_colours_in-1)*t2'),:); 

colormap(colour_map_transformed);

%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(9); clf;
fig.Name = 'PTC Scans: Intensity';
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

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = pad_data(data_hole_lt1, 0, 'lt1');
Z(end, end) = Z(end, end-1);
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = pad_data(data_hole_gt1, 0, 'gt1');
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: Before hole
[X, Y, Z] = pad_data(data_before_hole, 0, 'none');
surf(ax, X, Y, Z, MeshStyle='column');

% Surface: After hole
[X, Y, Z] = pad_data(data_after_hole, 0, 'none');
surf(ax, X, Y, Z, MeshStyle='column');

%-----------------------------------------%
%     Plot: Surface (One Level Lower)     %
%-----------------------------------------%
% % Surface: Hole (theta_old < 1)
% [X, Y, Z] = pad_data(data_hole_lt1, 1, 'lt1');
% Z(end, end) = Z(end, end-1);
% surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');
% 
% % Surface: Hole (theta_old > 1)
% [X, Y, Z] = pad_data(data_hole_gt1, 1, 'gt1');
% surf(ax, X, Y, Z, MeshStyle='column');
% 
% % Surface: Before hole
% [X, Y, Z] = pad_data(data_before_hole, 1, 'none');
% surf(ax, X, Y, Z, MeshStyle='column');
% 
% % Surface: After hole
% [X, Y, Z] = pad_data(data_after_hole, 1, 'none');
% surf(ax, X, Y, Z, MeshStyle='column');

%---------------------------%
%     Surface: Settings     %
%---------------------------%
% Shading of surface
shading(ax, 'interp');

% % Colorbar
% cbar = colorbar();

hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 25.0];
ax.ZAxis.Limits = [0.25, 3.0];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 5.0 : 25.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 2.5 : 25.0;

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
view(315, 15);
filename_out = '../fig8_I_PTC_surface_1.png';

% exportgraphics(fig, filename_out, ContentType='image', Resolution=750);

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
    [TO_min_val, max_idx] = min(min_theta_old);

    % Interpolate data
    theta_old_max = theta_old_read{max_idx};
  end

  if strcmp(theta_old_lt1_gt1, 'gt1')
    % Find theta_old_read with max theta_old value, interpolate so is the
    % same length as theta_old_max
    [TO_max_val, max_idx] = max(max_theta_old);

    % Interpolate data
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

      % Find all values less than this
      lt1_idx = theta_old_interp < min_val;

      % Set to the value
      theta_old_interp(lt1_idx) = min_val;

      % Set all NaN values to the first theta_new value
      theta_new_NaN = theta_new_read{i}(min_idx);

      % Find NaNs
      is_it_a_nan = isnan(theta_new_interp);

      % Find all NaNs below the minimum theta_old value
      theta_new_interp(is_it_a_nan(lt1_idx)) = theta_new_NaN;
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
      is_it_a_nan = isnan(theta_new_interp);

      % Find all NaNs above the maximum theta_old value
      theta_new_interp(is_it_a_nan) = theta_new_NaN;
      % theta_new_interp(is_it_a_nan(gt1_idx)) = theta_new_NaN;
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