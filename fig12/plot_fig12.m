% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig12_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Default colour order
colours = colororder();

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

% Setup figure
fig = figure(5); clf;
fig.Name = 'DTCs';
ax = gca();

% Axis dimensions
width = 6;
height = 8.5;

% Set figure size
set_figure_dimensions(width, height, scale=2);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%---------------------%
%     Plot: Patch     %
%---------------------%
% Fundamental domain
patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
    FaceAlpha=0.2, EdgeColor='none', ...
    HandleVisibility='off');

%--------------------%
%     Plot: DTCs     %
%--------------------%
% Plot linewidth
lw = 1.5;
% Plot colours
DTC_colours = {colours(9, :); colours(4, :); colours(5, :)};

for idx_A = 1 : 3
  % Read data
  x_plot = theta_perturb{idx_A};
  y_plot = theta_new{idx_A};

  % Plot DTC
  for offset = -1 : 1
    plot(ax, x_plot, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
    plot(ax, x_plot+1, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
  end
end
%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 0.5;

ax.YAxis.TickValues = -0.5 : 0.5 : 1.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.25 : 0.25 : 1.25;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\varphi_{\mathrm{d}}$');
% ylabel(ax, '$\vartheta_{\mathrm{n}}$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1];
ax.YAxis.Limits = [-0.25, 1.25];


%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../pdf/fig12_DTCs.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
