clear all; close all; clc;

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
fig = figure(7); clf;
fig.Name = 'Phase Reset in time: Intensity';
ax = gca();

% Axis dimensions
width = 3.8;
height = 1.5;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%------------------------------%
%     Plot: Phase Portrait     %
%------------------------------%
% Hold axes
hold(ax, 'on');

% Plot unerperturbed orbit
plot(ax, tbp_PO_plot, xbp_PO_plot, Color=[colours(3, :)], LineWidth=1.0);

% Plot segment 4
plot(ax, tbp4_run2, xbp4_run2(:, 3), Color=[0.0, 0.0, 0.0, 0.5], LineWidth=1.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
% ax.XAxis.Limits = [29, 30];
ax.XAxis.Limits = [27.5, 29.3];
ax.YAxis.Limits = [-0.1, 20];

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
ax.XAxis.TickValues = ax.XAxis.Limits;
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1) : 0.3 : ax.XAxis.Limits(2);

% Y-Axis
ax.YAxis.TickValues = 0.0 : 10 : 20.0;
ax.YAxis.MinorTickValues = 0.0 : 5 : 20.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$t / T_{\Gamma}$');
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
filename_out = '../pdf/fig4b3_zoom.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
