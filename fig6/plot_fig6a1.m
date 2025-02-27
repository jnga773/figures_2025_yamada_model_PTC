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
% Plot colours
plot_colours = {'#bcbd22';    % Green-Yellow
                '#e66119';    % Orange
                '#8c564b';    % Brown
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5';    % Purple
                '#1f9ece'};   % Cyan

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
width = 7.5;
height = 4.0;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%--------------%
%     Plot     %
%--------------%
% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), Color=colours(3, :), ...
      LineWidth=2.0);

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), Marker='o', MarkerSize=7.5, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%----------------------------%
%     Plot: Perturbation     $
%----------------------------%
lw = 1.0;

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % Plot
  plot3(ax, smooth(xbp_PO(:, 1)+A_perturb(idx)), smooth(xbp_PO(:, 2)), smooth(xbp_PO(:, 3)), ...
        Color=plot_colours{idx}, LineStyle='-', LineWidth=lw)%, ...
        % Marker='none', MarkerFaceColor=colour(1:3), MarkerEdgeColor=colour(1:3));
end

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 7.0];
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, 20];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 2.0 : 8.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 7.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 2.0 : 4.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 1.0 : 4.0;

% Z-Axis
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickValues = 0.0 : 5.0 : 20.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 0.0 : 2.5 : 20.0;

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

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% 3D plot view
% view(45, 10.0);
view(45, 6.0);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig6a1_phase_portrait.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
