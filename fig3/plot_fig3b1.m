clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read from data structure file
load('../data_files/fig3_data.mat');

%-------------------------%
%     Read Parameters     %
%-------------------------%
% Print parameters to console
fprintf('A_perturb(1) = %.4f\n', A_perturb_run1);
fprintf('A_perturb(2) = %.4f\n', A_perturb_run2);
fprintf('theta_old = %.4f\n', theta_old);

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(5); clf;
fig.Name = 'Phase Reset Phase Portrait (2D)';
ax = gca();

% Axis dimensions
width = 3.7;
height = 3.7;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%------------------------------%
%     Plot: Phase Portrait     %
%------------------------------%
% Hold axes
hold(ax, 'on');

% Plot original periodic orbit
plot(ax, xbp_PO(:, 1), xbp_PO(:, 3), ...
     Color=colours(3, :), ...
     LineWidth=2.0);

% Plot segment 4
plot(ax, xbp4_run2(:, 1), xbp4_run2(:, 3), ...
     Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=0.5);

% Plot equilibrium point
plot(ax, xpos(1), xpos(3), ...
     Marker='o', MarkerSize=4.0, ...
     MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% Plot theta_old point on gamma
plot(ax, xbp_PO(1, 1), xbp_PO(1, 3), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor=colours(3, :), MarkerEdgeColor='k');
     
% Plot start point of segment 4
plot(ax, xbp4_run2(1, 1), xbp4_run2(1, 3), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 7];
ax.YAxis.Limits = [-0.5, 50];

%------------------------------%
%     Axis Ticks: Settings     %
%------------------------------%
ax.XAxis.TickDirection = 'in';
ax.XAxis.MinorTick = 'on';
ax.YAxis.TickDirection = 'in';
ax.YAxis.MinorTick = 'on';

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 2 : 8;
ax.XAxis.MinorTickValues = 0.0 : 1 : 8.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 10 : 50.0;
ax.YAxis.MinorTickValues = 0.0 : 5 : 50.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$G$');
% ylabel(ax, '$I$');

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
filename_out = '../pdf/fig3b1_portrait2.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
