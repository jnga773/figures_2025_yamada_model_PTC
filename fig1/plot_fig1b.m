clear all; close all; clc;

%-------------------------------------------------------------------------%
%%                         Read Bifurcation Data                         %%
%-------------------------------------------------------------------------%
load('../data_files/fig1_data.mat');

%%
%-----------------------------------------------------------------------%
%                          Plot: Big Picture                            %
%-----------------------------------------------------------------------%
% Plot colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Yamada - Bifurcations';
ax = gca();

% Axis dimensions
% width = 7.5;
% height = 4;
width = 4.5;
height = 2.5;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%--------------%
%     Plot     %
%--------------%
% Hold axis
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
     Color=colours(4, :));

% % Plot: Neutral saddle-Node
% plot(ax, A_NSA, gamma_NSA, ...
%      LineStyle='--', LineWidth=1.5, ...
%      Color=colours(4, :));

% Plot: Homoclinics (approximate)
plot(ax, A_L_approx, gamma_L_approx, ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colours(1, :));

% % Plot: Homoclinics (Lin's method)
% plot(ax, A_L_lins, gamma_L_lins, ...
%      LineStyle='-', LineWidth=1.5, ...
%      Color=colours(1, :));

% Plot: Double limit cycles
plot(ax, A_D, gamma_D, ...
     LineStyle='-', LineWidth=1.5, ...
     Color=colours(3, :));

% Add dot for phase resetting parameters
plot(ax, 7.4, 3.5e-2, Marker='pentagram', MarkerFaceColor='k', MarkerEdgeColor='k');

% Turn off axis hold
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [6.65, 7.5];
ax.YAxis.Limits = [0.03, 0.08];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues =[];

% Y-Axis
ax.YAxis.TickValues = [];

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
filename_out = '../pdf/fig1b_zoomed_bifurcation_diagram.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
