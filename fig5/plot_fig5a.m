% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
% Read data from .mat file
load('../data_files/fig5_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
plot_colours = {'#eb5e00';    % Orange
                '#bf42f5'};   % Purple

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
height = 7.0;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%----------------------------%
%     Plot: Other Things     %
%----------------------------%
% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), ...
      FaceAlpha=0.2, EdgeColor='none');

% Plot diagonal lines
plot(ax, [0, 1], [0, 1], LineStyle='-', LineWidth=1.5, ...
     Color=colours(3, :));
plot(ax, [0, 1], [-1, 0], LineStyle='-', LineWidth=1.5, ...
     Color=colours(3, :));
plot(ax, [0, 1], [1, 2], LineStyle='-', LineWidth=1.5, ...
     Color=colours(3, :));

% Grey lines at theta_old = 0 and 0.3
xline(ax, 0.3, LineStyle='-', LineWidth=1, ...
      Color=[0, 0, 0, 0.5]);

%-------------------%
%     Plot: PTC     %
%-------------------%
% Plot PTC: run1
for i = 1 : length(theta_old)
  plot(ax, theta_old{i}, theta_new{i}-1, ...
       LineStyle='-', LineWidth=1.5, ...
       Color=plot_colours{i});
  plot(ax, theta_old{i}, theta_new{i}, ...
       LineStyle='-', LineWidth=1.5, ...
       Color=plot_colours{i});
  plot(ax, theta_old{i}, theta_new{i}+1, ...
       LineStyle='-', LineWidth=1.5, ...
       Color=plot_colours{i});
end

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
ax.YAxis.Limits = [-0.1, 1.25];

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
filename_out = '../pdf/fig5_G_two_PTCs.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
