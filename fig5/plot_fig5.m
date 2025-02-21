% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
% Read data from .mat file
load('../data_files/fig5_data.mat');

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
plot_colours = {[250, 250,  42] ./ 255;
                [148, 103, 189] ./ 255};

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
% xline(ax, theta_old_run1, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);
xline(ax, 0.3, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);

%-------------------%
%     Plot: PTC     %
%-------------------%
% Plot PTC: run1
for i = 1 : length(theta_old)
  plot(ax, theta_old{i}, theta_new{i}-1, Color=plot_colours{i}, LineStyle='-', LineWidth=1.5);
  plot(ax, theta_old{i}, theta_new{i}, Color=plot_colours{i}, LineStyle='-', LineWidth=1.5);
  plot(ax, theta_old{i}, theta_new{i}+1, Color=plot_colours{i}, LineStyle='-', LineWidth=1.5);
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
filename_out = '../pdf/fig5_G_PTC_single.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
