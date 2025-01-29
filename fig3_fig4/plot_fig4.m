clear all; close all; clc;

%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
% Load data
load('../plot_mat_files/fig4_data.mat');

%%
%-------------------------------------------------------------------------%
%                       Plot Phase Transition Curve                       %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Single PTC';

% Figure dimensions
fig.Units = 'centimeters';
fig.Position = [5, 5, 6.0, 6.0];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

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
xline(ax, 0, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);
xline(ax, 0.3, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);

%-------------------%
%     Plot: PTC     %
%-------------------%
% Plot PTC: run1
plot(ax, theta_old_run1, theta_new_run1, Color='k', LineStyle='-', LineWidth=1.5);
plot(ax, theta_old_run1, theta_new_run1+1, Color='k', LineStyle='-', LineWidth=1.5);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

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

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% % Turn off all axis labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\\theta_{\\mathrm{o}}$';
ax.YAxis.Label.String = '$\\theta_{\\mathrm{n}}$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../images/fig4_G_PTC_single.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
