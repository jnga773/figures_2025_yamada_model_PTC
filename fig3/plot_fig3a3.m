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
fprintf('theta_old    = %.4f\n', theta_old);
fprintf('A_perturb(1) = %.4f\n', A_perturb_run1);
fprintf('theta_new(1) = %.4f\n', theta_new_run1);
fprintf('A_perturb(2) = %.4f\n', A_perturb_run2);
fprintf('theta_new(2) = %.4f\n', theta_new_run2);

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO  = '#2ca02c';
% Perturbed orbit colour
colour_PR  = [0.0, 0.0, 0.0, 0.5];

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(4); clf;
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
plot(ax, tbp_PO_plot, xbp_PO_plot, ...
     Color=colour_PO, LineWidth=1.0);

% Plot segment 4
plot(ax, tbp4_run1, xbp4_run1(:, 3), ...
     Color=colour_PR, LineWidth=1.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [29, 29.5];
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
ax.XAxis.MinorTickValues = ax.XAxis.Limits(1) : 0.1 : ax.XAxis.Limits(2);

% Y-Axis
ax.YAxis.TickValues = ax.YAxis.Limits;
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
filename_out = './fig3a3_zoom.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
