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
intersection.theta_old = 0.3176;
intersection.A_perturb = 0.5576;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
plot_colours = {'#bcbd22';    % Green-Yellow
                '#e66119';    % Orange
                '#8c564b';    % Brown
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5';    % Purple
                '#1f9ece'};   % Cyan

%-----------------------------------%
%     Sort Out Single Plot Data     %
%-----------------------------------%
plot_A_perturb = [0.05, 0.1, 0.15, 0.55, 1.0, 1.5, 2.0];

% Find plotting indices
plot_idx = zeros(length(plot_A_perturb), 1);
for i = 1 : length(plot_A_perturb)
  plot_idx(i) = find(round(A_perturb, 3) == plot_A_perturb(i));
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

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(1); clf;
fig.Name = 'PTC Scans: Gain';
ax = gca();

% Axis dimensions
width = 8.0;
height = 6.4;

% Set figure size
set_figure_dimensions(width, height, scale=4);

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
colormap(scale_colour_map(2.0));

% Shading of surface
shading(ax, 'interp');

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

%--------------------------%
%     Plot: PTC Curves     %
%--------------------------%
% Linewidth
lw = 3.0;

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);
  plot3(ax, theta_old_plot{i}, A_perturb_plot{i}, theta_new_plot{i}, ...
        LineWidth=lw, LineStyle='-', ...
        Color=plot_colours{i});
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
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};
% ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% axis(ax, 'off');
% set(gca, 'Color', 'none');

%---------------------%
%     Save Figure     %
%---------------------%
view(315, 15);

filename_out = '../pdf/fig7a_G_PTC_surface_1.png';
% exportgraphics(fig, filename_out, ContentType='image', Resolution=1000);

% filename_out = '../pdf/fig7a_G_PTC_surface_1.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
