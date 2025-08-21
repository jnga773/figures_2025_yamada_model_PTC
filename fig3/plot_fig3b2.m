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
fprintf('A_perturb(1) = %.4f\n', A_perturb_run1);
fprintf('A_perturb(2) = %.4f\n', A_perturb_run2);
fprintf('theta_old = %.4f\n', theta_old);

%----------------------%
%     Plot Colours     %
%----------------------%
% Matplotlib colours
colours = {'#1f77b4';  % blue
           '#ff7f0e';  % orange
           '#2ca02c';  % green
           '#d62728';  % red
           '#9467bd';  % purple
           '#8c564b';  % brown
           '#e377c2';  % pink
           '#7f7f7f';  % gray
           '#bcbd22';  % yellow-green
           '#17becf'   % cyan
           };

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(6); clf;
fig.Name = 'Phase Reset in time: Intensity';
ax = gca();

% Axis dimensions
width = 7.8;
height = width / 3;

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
     Color=[colours{3}], ...
     LineWidth=1.0);

% Plot segment 4
max_idx = max(find(tbp4_run2 < 12.0));
plot(ax, tbp4_run2(1:max_idx), xbp4_run2(1:max_idx, 3), ...
     Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=1.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-0.2, 12];
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
ax.XAxis.TickValues = 0.0 : 2.0 : 12.0;
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 12.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 10 : 50.0;
ax.YAxis.MinorTickValues = 0.0 : 5 : 50.0;

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
filename_out = './fig3b2_time2.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
