% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig2_data.mat', 'Wq_s');
load('../data_files/fig6_data.mat');

% List of perturbations to plot
plot_idx = 1:4;
% plot_idx = 4:7;

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO    = '#2ca02c';
% Stable manifold colour
colour_Wsq   = '#1f77b4';

% Plot colours
plot_colours = {'#92b700';    % Green-Yellow
                '#e6b400';    % Yellow
                '#eb5e00';    % Orange
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5';    % Purple
                '#1f9ece'};   % Cyan

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(2); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
width = 3.0;
height = 1.6;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%--------------%
%     Plot     %
%--------------%
% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
      Color=colour_PO, LineWidth=2.0);

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours_Wsq, LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), ...
      Marker='o', MarkerSize=4, ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%----------------------------%
%     Plot: Perturbation     %
%----------------------------%
lw = 2.0;

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % Plot
  plot3(ax, smooth(xbp_PO(:, 1)+A_perturb(idx)), smooth(xbp_PO(:, 2)), smooth(xbp_PO(:, 3)), ...
        Color=plot_colours{idx}, LineStyle='-', LineWidth=lw);
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.75, 2.5];
ax.YAxis.Limits = [0.0, 1.5];
ax.ZAxis.Limits = [0.0, 5];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = [];
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = [];

% Y-Axis
ax.YAxis.TickValues = [];
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = [];

% Z-Axis
ax.ZAxis.TickValues = [];
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = [];

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$G$');
% ylabel(ax, '$Q$');
% zlabel(ax, '$I$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};
ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% 3D plot view
% view(45, 10.0);
view(-30, 6.0);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = './fig6a1_inset.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
