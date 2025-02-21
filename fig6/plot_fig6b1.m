% close all; clear all; clc

%-------------------------------------------------------------------------%
%%                               Read Data                               %%
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig6_data.mat');

% Data indices to plot
plot_idx = 1:4;
% plot_idx = 4:7;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
plot_colours = {'#bcbd22';    % Green-Yellow
                '#d8a400';    % Yellow-Orange
                '#e66119';    % Orange
                '#d62728';    % Red
                '#a12b6f';    % Purple-Magenta
                '#5a5fc8';    % Blue-Violet
                '#1f9ece'};   % Cyan

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

% Setup figure
fig = figure(5); clf;
fig.Name = 'PTCs';
ax = gca();

% Axis dimensions
width = 4.5;
height = 9.0;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%----------------------------%
%     Plot: Other Things     %
%----------------------------%
% Vertical line at intersection point
xline(ax, 0.3176, Color=[0.0 0.0 0.0], LineStyle='-', LineWidth=1.0);

% Plot diagonal line% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='-', Color=[ 44, 160,  44] ./ 255, LineWidth=1.5, ...
     HandleVisibility='off');

% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

%--------------------%
%     Plot: PTCs     %
%--------------------%
% Linewidth
lw = 1.5;

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % Plot
  plot(ax, theta_old{idx}, theta_new{idx}, Color=plot_colours{idx}, LineStyle='-');
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
ax.DataAspectRatio = [1, 1, 1];

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = -0.5 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -0.5 : 0.25 : 1.0;

ax.YAxis.TickValues = -0.5 : 0.5 : 2.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.5 : 0.25 : 2.5;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\theta_{\mathrm{o}}$');
% ylabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1.0];
ax.YAxis.Limits = [-0.25, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../pdffig6b1_G_PTCs.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
