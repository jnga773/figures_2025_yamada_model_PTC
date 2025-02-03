clear all; close all; clc;

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
fig = figure(1); clf;
fig.Name = 'Single PTC';

% Figure dimensions
% fig.Units = 'centimeters';
fig.Units = 'inches';
fig.Position = [5, 5, 6.0, 6.0];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% % Axis setup: Manual padding
% ax = gca();
% ax.Position = [0.01, 0.01, 0.98, 0.98];

% Axis setup: Tiled layout
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

% Set fontsizes
ax.FontSize = 9;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%----------------------------%
%     Plot: Other Things     %
%----------------------------%
% Straight line for intersection point theta_old
% xline(ax, 0.3176, Color=[0.0 0.0 0.0], LineStyle='-', LineWidth=2.5, ...
%       HandleVisibility='off')

% % Points for theta_old
% idx_3 = find(round(theta_old_run1, 4) == 0.3);
% plot(ax, [theta_old_run1(1), theta_old_run1(idx_3)], [theta_new_run1(1), theta_new_run1(idx_3)], ...
%      LineStyle='none', Marker='o', MarkerFaceColor='r', MarkerEdgeColor='r', MarkerSize=5);

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
xlabel(ax, '$\theta_{\mathrm{o}}$');
ylabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all tick labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
% filename_out = '../fig4.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
