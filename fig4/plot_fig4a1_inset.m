% clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read from data structure file
load('../data_files/fig4_data.mat');

%-------------------------%
%     Read Parameters     %
%-------------------------%
% Print parameters to console
fprintf('theta_old    = %.4f\n', theta_old);
fprintf('A_perturb(1) = %.4f\n', A_perturb_run1);
fprintf('theta_new(1) = %.4f\n', theta_new_run1);
fprintf('A_perturb(2) = %.4f\n', A_perturb_run2);
fprintf('theta_new(2) = %.4f\n', theta_new_run2);

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO = '#2ca02c';
% Perturbed orbit colour
colour_PR  = [0.0, 0.0, 0.0, 0.5];

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(2); clf;
fig.Name = 'Phase Reset Phase Portrait (2D)';
ax = gca();

% Axis dimensions
width = 1.5;
height = 1.0;

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
     Color=colour_PO, LineWidth=2.0);

% Plot segment 4
plot(ax, xbp4_run1(:, 1), xbp4_run1(:, 3), ...
     Color=colour_PR, LineWidth=0.5);

% Plot equilibrium point
plot(ax, xpos(1), xpos(3), ...
     Marker='o', MarkerSize=4.0, ...
     MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% Plot theta_old point on gamma
plot(ax, xbp3_run1(1, 1), xbp3_run1(1, 3), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor=colour_PO, MarkerEdgeColor='k');

% Plot start point of segment 4
plot(ax, xbp4_run1(1, 1), xbp4_run1(1, 3), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.96, 1.15];
ax.YAxis.Limits = [2.9, 3.3];

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
ax.XAxis.TickValues = [];
ax.XAxis.MinorTickValues = [];

% Y-Axis
ax.YAxis.TickValues = [];
ax.YAxis.MinorTickValues = [];

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
filename_out = './fig4a1_inset.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
