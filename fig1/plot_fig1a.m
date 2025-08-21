% clear all; close all; clc;

%-------------------------------------------------------------------------%
%%                         Read Bifurcation Data                         %%
%-------------------------------------------------------------------------%
load('../data_files/fig1_data.mat');

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
%-----------------------------------------------------------------------%
%                          Plot: Big Picture                            %
%-----------------------------------------------------------------------%
% Setup figure
fig = figure(1); clf;
fig.Name = 'Yamada - Bifurcations';
ax = gca();

% Axis dimensions
width = 7.5;
height = 4;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%--------------%
%     Plot     %
%--------------%
hold(ax, 'on');

% Plot: Saddle-node bifurcations
plot(ax, A_S, gamma_S, ...
     LineStyle='-', LineWidth=1.5, ...
     Color='k');

% Plot: Transcritical bifurcations
plot(ax, A_T, gamma_T, ...
     LineStyle='-', LineWidth=1.5, ...
     Color='k');

% Plot: Hopf birfurcations
plot(ax, A_H, gamma_H, ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colours{4});

% % Plot: Neutral saddle-Node
% plot(ax, A_NSA, gamma_NSA, ...
%      LineStyle='--', LineWidth=1.5, ...
%      Color=colours{4});

% Plot: Homoclinics (approximate)
plot(ax, A_L_approx, gamma_L_approx, ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colours{1});

% % Plot: Homoclinics (Lin's method)
% plot(ax, A_L_lins, gamma_L_lins, ...
%      LineStyle='-', LineWidth=1.5, ...
%      Color=colours{1});

% Plot: Double limit cycles
plot(ax, A_D, gamma_D, ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colours{3});

% % Add dot for phase resetting parameters
% plot(ax, 7.4, 3.5e-2, Marker='pentagram', MarkerFaceColor='k', MarkerEdgeColor='k');

% % Plot box for zoom
% rectangle(ax, Position=[6.65, 0.03, 7.5-6.65, 0.08-0.03], LineWidth=0.8, ...
%           EdgeColor='k');
% % ax.XAxis.Limits = [6.65, 7.5];
% % ax.YAxis.Limits = [0.03, 0.08];

% Turn off axis hold
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [5.8, 11];
ax.YAxis.Limits = [0.0, 0.25];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 6 : 1 : 11;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 6 : 0.5 : 11;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 0.1 : 0.30;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.00 : 0.05 : 0.30;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$A_{\mathrm{p}}$');
% ylabel(ax, '$\gamma$');

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
filename_out = './fig1a_bifurcation_diagram.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
