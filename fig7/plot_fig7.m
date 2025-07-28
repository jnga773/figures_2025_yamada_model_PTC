clear all; close all; clc;

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
intersection.theta_old = 0.32601;
intersection.A_perturb = 0.54353;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
plot_colours = {'#92b700';    % Green-Yellow
                '#e6b400';    % Yellow
                '#eb5e00';    % Orange
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5';    % Purple
                '#1f9ece'};   % Cyan

%-----------------------------------%
%     Sort Out Single Plot Data     %
%-----------------------------------%
plot_A_perturb = [0.05, 0.1, 0.15, 0.5427, 1.0, 1.5, 2.0];

% Find plotting indices
plot_idx = zeros(length(plot_A_perturb), 1);
for i = 1 : length(plot_A_perturb)
  plot_idx(i) = find(round(A_perturb, 4) == plot_A_perturb(i));
end

% Empty cells for plotting data
theta_old_plot = cell(1, length(plot_idx));
theta_new_plot = cell(1, length(plot_idx));
A_perturb_plot = cell(1, length(plot_idx));

for i = 1 : length(plot_A_perturb)
  % Plot data index
  idx = plot_idx(i);

  % Logical checks
  lt1_check = false;
  gt1_check = false;

  if round(theta_old_lt1{idx}(end) - theta_old_lt1{idx}(1), 3) == 1.0
    lt1_check = true;
  end
  if round(theta_old_gt1{idx}(end) - theta_old_gt1{idx}(1), 3) == 1.0
    gt1_check = true;
  end

  % Grab data
  if lt1_check && gt1_check
    % Both checks are true; take the lt1 data
    theta_old_plot{i} = theta_old_lt1{idx};
    theta_new_plot{i} = theta_new_lt1{idx};
    A_perturb_plot{i} = A_perturb(idx) * ones(1, length(theta_old_lt1{idx}));
  elseif lt1_check && ~gt1_check
    % Only the lt1 check is true
    theta_old_plot{i} = theta_old_lt1{idx};
    theta_new_plot{i} = theta_new_lt1{idx};
    A_perturb_plot{i} = A_perturb(idx) * ones(1, length(theta_old_lt1{idx}));
  elseif ~lt1_check && gt1_check
    % Onle the gt1 check is true
    theta_old_plot{i} = theta_old_gt1{idx}-1.0;
    theta_new_plot{i} = theta_new_gt1{idx};
    A_perturb_plot{i} = A_perturb(idx) * ones(1, length(theta_old_gt1{idx}));
  else
    % Neither of the checks are true
    theta_old_plot{i} = [theta_old_gt1{idx}-1.0, nan, theta_old_lt1{idx}];
    theta_new_plot{i} = [theta_new_gt1{idx}, nan, theta_new_lt1{idx}];
    A_perturb_plot{i} = A_perturb(idx) * ones(1, length([theta_old_gt1{idx}, nan, theta_old_lt1{idx}]));
  end
end

% Shift data_hole_lt1 start data point
data_hole_lt1 = fix_gap(data_before_hole, data_hole_lt1);

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(1); clf;
fig.Name = 'PTC Scans: Gain';
ax  = axes;

% Axis dimensions
width = 7.8;
height = 6.0;

% Set figure size
set_figure_dimensions(width, height, scale=1);

% Set axis linewidth
ax.LineWidth = 0.8;

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
%     Surface Settings     %
%--------------------------%
% Set colour map
cmap = colormap(scale_colour_map(2.0));

% Shading of surface
shading(ax, 'interp');

% % Lighting
% light(ax, Style='infinite', Position=[intersection.theta_old, intersection.A_perturb, 20]);

% Facealpha
facealpha = 0.75;

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = pad_data(data_hole_lt1, 0, 'lt1');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z+1.0, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = pad_data(data_hole_gt1, 0, 'gt1');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

% Surface: Before hole
[X, Y, Z] = pad_data(data_before_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z-1, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

% Surface: After hole
[X, Y, Z] = pad_data(data_after_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z+1, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

%--------------------------%
%     Plot: PTC Curves     %
%--------------------------%
% Linewidth
lw = 1.0;
% 
% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);
  plot3(ax, theta_old_plot{i}, A_perturb_plot{i}, theta_new_plot{i}, ...
        LineWidth=lw, LineStyle='-', ...
        Color=plot_colours{i});
  if i <= 3
    plot3(ax, theta_old_plot{i}, A_perturb_plot{i}, theta_new_plot{i}-1, ...
          LineWidth=lw, LineStyle='-', ...
          Color=plot_colours{i});
  elseif i > 3
    plot3(ax, theta_old_plot{i}, A_perturb_plot{i}, theta_new_plot{i}+1, ...
          LineWidth=lw, LineStyle='-', ...
          Color=plot_colours{i});
  end
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 2.0];
ax.ZAxis.Limits = [0.0, 3.25];

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

% axis(ax, 'off');

%---------------------%
%     Save Figure     %
%---------------------%
view(315, 15);

% filename_out = '../pdf/fig7_G_PTC_surface.png';
% exportgraphics(fig, filename_out, ContentType='image', Resolution=1000);

% filename_out = '../pdf/fig7_G_PTC_surface.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');

%-------------------------------------------------------------------------%
%%                               FUNCTION                                %%
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

function data_out = fix_gap(data_before_hole, data_hole_lt1)
  % Fills in the gap between data_hole_lt1 and the level below of
  % data_before_hole.

  % Get data from end of data_before_hole
  A_perturb_before_hole = data_before_hole.A_perturb(end);
  theta_old_before_hole = data_before_hole.theta_old{end};
  theta_new_before_hole = data_before_hole.theta_new{end} - 1.0;

  % Get data from start of data_hole_lt1
  A_perturb_hole = data_hole_lt1.A_perturb(1);
  theta_old_hole = data_hole_lt1.theta_old{1};
  theta_new_hole = data_hole_lt1.theta_new{1};

  % Find nearest point in theta_old_before_hole to the start point of
  % theta_old_hole
  [min_val, min_idx] = min(abs(theta_old_before_hole - theta_old_hole(1)));

  % Append the data
  A_perturb_out = [A_perturb_before_hole, data_hole_lt1.A_perturb(2:end)];
  theta_old_out = {theta_old_before_hole(min_idx:end), data_hole_lt1.theta_old{2:end}};
  theta_new_out = {theta_new_before_hole(min_idx:end), data_hole_lt1.theta_new{2:end}};

  data_out.A_perturb = A_perturb_out;
  data_out.theta_old = theta_old_out;
  data_out.theta_new = theta_new_out;

end

function colourmap_out = scale_colour_map(scale_factor)
  % colourmap_out = scale_colour_map(scale_factor)
  % 
  % Rescale the colour map.

  n_colours_in = 2048;
  n_colours_out = 1024;

  % Get colour map
  colour_map = parula(n_colours_in);

  % Create a linear ramp the size of the colormap we actually want
  t = linspace(0,1,n_colours_out);
  % Apply whatever transform you like to the ramp
  t2 = t .^ scale_factor;

  % Use that to scale the big linear colormap into the small stretched one.
  colourmap_out = colour_map(1+floor((n_colours_in-1)*t2'),:);
end
