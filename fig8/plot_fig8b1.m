% clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read from data structure file
load('../data_files/fig8_data.mat');

%-------------------------%
%     Read Parameters     %
%-------------------------%
% Print parameters to console
fprintf('A_perturb(1) = %.4f\n\n', A_perturb_run1);
fprintf('A_perturb(2) = %.4f\n\n', A_perturb_run2);
fprintf('theta_old = %.4f\n\n', theta_old);

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(3); clf;
fig.Name = 'Phase Reset in time: Intensity';
ax = gca();

% Axis dimensions
width = 7.8;
height = width / 3;

% Add set_figure_dimensions() function to path
% addpath('../');

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
max_idx = max(find(tbp_PO_plot < 12.0));
plot(ax, tbp_PO_plot(1:max_idx), xbp_PO_plot(1:max_idx), ...
     Color=[colours(3, :)], ...
     LineWidth=1.0);

% Plot segment 4
max_idx = max(find(tbp4_run1 < 12.0));
plot(ax, tbp4_run1(1:max_idx), xbp4_run1(1:max_idx, 3), ...
     Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=1.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-0.2, 12];
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
ax.XAxis.TickValues = 0.0 : 2.0 : 12.0;
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 12.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 5 : 20.0;
ax.YAxis.MinorTickValues = 0.0 : 2.5 : 20.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$t / T_{\Gamma}$');
% ylabel(ax, '$I$');

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
% filename_out = '../fig8b1_time1.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
