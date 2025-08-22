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

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
colour_G = '#1f77b4';
colour_Q = '#ff7f0e';
colour_I = '#2ca02c';

%%
%-------------------------------------------------------------------------%
%                          Plot: Temporal Trace                           %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(1); clf;
fig.Name = 'Temporal Trace of Periodic Orbit';
ax = gca();

% Axis dimensions
width = 7.5;
height = 3.0;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot: Gain
plot(ax, tbp_plot, xbp_plot(:, 1), ...
     LineStyle=':', LineWidth=1.5, ...
     Color=colour_G);

% Plot: Absorption
plot(ax, tbp_plot, xbp_plot(:, 2), ...
     LineStyle='--', LineWidth=1.5, ...
     Color=colour_Q);

% Plot: Intensity
plot(ax, tbp_plot, xbp_plot(:, 3), ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colour_I);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 2];
ax.YAxis.Limits = [0.0, 16.5];

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
ax.YAxis.TickValues = 0.0 : 4.0 : 16.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 2.0 : 16.0;

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
filename_out = './fig2a_periodic_orbit_temporal_trace.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
