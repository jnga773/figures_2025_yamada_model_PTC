% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat');

%------------------------%
%    Add two periods     %
%------------------------%
% Period of periodic orbit
tbp = tbp_PO / T_PO;

% Copy original periodic orbit data
xbp_plot = [xbp_PO(1:end-1, :); xbp_PO];
tbp_plot = [tbp(1:end-1); 1 + tbp];

%%
%-------------------------------------------------------------------------%
%                          Plot: Temporal Trace                           %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(2); clf;
fig.Name = 'Temporal Trace of Periodic Orbit';
ax = gca();

% Axis dimensions
width = 7.5;
height = 3.0;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot original periodic orbit
plot(ax, tbp_plot, xbp_plot(:, 1), Color=colours(1, :), LineStyle=':', LineWidth=1.5);
plot(ax, tbp_plot, xbp_plot(:, 2), Color=colours(2, :), LineStyle='--', LineWidth=1.5);
plot(ax, tbp_plot, xbp_plot(:, 3), Color=colours(3, :), LineStyle='-', LineWidth=1.5);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 2];
ax.YAxis.Limits = [-0.05, 20];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 0.5 : 2.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 2.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 5.0 : 25.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 2.5 : 25.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$t / T_{\Gamma}$');
% ylabel(ax, '$G/Q/I$');

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
filename_out = '../pdf/fig2b_periodic_orbit_temporal_trace.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
