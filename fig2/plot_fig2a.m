clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
%----------------------------------%
%     Read Data from .mat File     %
%----------------------------------%
load('../data_files/fig2_data.mat');

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';

% Figure dimensions
% fig.Units = 'centimeters';
fig.Units = 'inches';
fig.Position = [5, 5, 6, 6];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), Color=colours(3, :), ...
      LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), Marker='diamond', MarkerSize=10, ...
      MarkerFaceColor='r', MarkerEdgecolor='r');

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), LineWidth=2.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, 21];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 2.0 : 6.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 6.0;

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
ax.XAxis.Label.String = '$G$';
ax.YAxis.Label.String = '$Q$';
ax.ZAxis.Label.String = '$I$';

% Turn off all axis labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};
% ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% 3D plot view
% view(45, 10);
view(45, 6);

%---------------------%
%     Save Figure     %
%---------------------%
% filename_out = '../fig2a_phase_portrait.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
