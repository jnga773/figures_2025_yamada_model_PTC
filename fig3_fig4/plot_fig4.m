% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
% Read data from .mat file
load('../data_files/fig4_data.mat');

%%
%-------------------------------------------------------------------------%
%                       Plot Phase Transition Curve                       %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(5); clf;
fig.Name = 'Single PTC';
ax = gca();

% Axis dimensions
width = 6.0;
height = 6.0;

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
% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='-', Color=colours(3, :), LineWidth=1.5, ...
     HandleVisibility='off');

% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

% Grey lines at theta_old = 0 and 0.3
xline(ax, theta_old_run1, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);
xline(ax, theta_old_run2, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);

%-------------------%
%     Plot: PTC     %
%-------------------%
% Plot PTC: run1
plot(ax, theta_old_plot, theta_new_plot, Color='k', LineStyle='-', LineWidth=1.5);
plot(ax, theta_old_plot, theta_new_plot+1, Color='k', LineStyle='-', LineWidth=1.5);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
ax.DataAspectRatio = [1, 1, 1];

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-0.0, 1.0];
ax.YAxis.Limits = [-0.1, 1.0];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = -1.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -1.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\theta_{\mathrm{o}}$');
% ylabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig4_G_PTC_single.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
