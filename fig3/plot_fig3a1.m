% clear all; close all; clc;

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
fprintf('A_perturb = %.4f\n\n', A_perturb);

fprintf('theta_old(1) = %.4f\n', theta_old_run1);
fprintf('theta_old(2) = %.4f\n\n', theta_old_run2);

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Phase Reset Phase Portrait (2D)';
ax = gca();

% Axis dimensions
width = 3.5;
height = 3.5;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

%------------------------------%
%     Plot: Phase Portrait     %
%------------------------------%
% Hold axes
hold(ax, 'on');

% Plot segment 4
plot(ax, xbp4_run1(:, 1), xbp4_run1(:, 3), Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=1.0, DisplayName='Segment 4');

% Plot original periodic orbit
plot(ax, xbp_PO(:, 1), xbp_PO(:, 3), Color=colours(3, :), ...
     LineWidth=2.0, DisplayName='$\Gamma$');

% Plot equilibrium point
plot(ax, xpos(1), xpos(3), Marker='o', MarkerSize=4.0, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 5];
ax.YAxis.Limits = [-0.1, 20];

%------------------------------%
%     Axis Ticks: Settings     %
%------------------------------%
ax.XAxis.TickDirection = 'in';
ax.XAxis.MinorTick = 'on';
ax.YAxis.TickDirection = 'in';
ax.YAxis.MinorTick = 'on';

%---------------------------------%
%     Axis Ticks: ax1 and ax2     %
%---------------------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 1 : 5.0;
ax.XAxis.MinorTickValues = 0.0 : 0.5 : 5.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 5 : 20.0;
ax.YAxis.MinorTickValues = 0.0 : 2.5 : 20.0;

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
filename_out = '../pdf/fig3a1_portrait1.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
