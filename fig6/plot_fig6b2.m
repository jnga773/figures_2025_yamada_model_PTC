close all; clear all; clc

%-------------------------------------------------------------------------%
%%                               Read Data                               %%
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig6_data.mat');

% Data indices to plot
% plot_idx = 1:4;
plot_idx = 4:7;

% Default line colours
colours = colororder();

% Plot colours
plot_colours = {colours(2, :), colours(4, :), colours(5, :), ...
                colours(6, :), ...
                colours(7, :), colours(8, :), colours(9, :)};

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(1); clf;
fig.Name = 'PTC Scan';

% Figure dimensions
fig.Units = 'centimeters';
fig.Position = [5, 5, 6, 12];
% fig.Position = [3, 3, 2, 5.5];
% fig.Position = [5, 5, 5, 10];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%----------------------------%
%     Plot: Other Things     %
%----------------------------%
% Vertical line at intersection point
xline(ax, 0.3176, Color=[0.0 0.0 0.0], LineStyle='-', LineWidth=1.0);

% Plot diagonal line% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='-', Color=colours(3, :), LineWidth=1.5, ...
     HandleVisibility='off');

% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

%--------------------%
%     Plot: PTCs     %
%--------------------%
% Linewidth
lw = 1.5;

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % Plot
  plot(ax, theta_old{idx}, theta_new{idx}, Color=plot_colours{idx}, LineStyle='-');
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = -0.5 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -0.5 : 0.25 : 1.0;

ax.YAxis.TickValues = -0.5 : 0.5 : 2.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.5 : 0.25 : 2.5;

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% Turn off all axis labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

xlabel(ax, '$\theta_{\mathrm{o}}$');
ylabel(ax, '$\theta_{\mathrm{n}}$');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1.0];
ax.YAxis.Limits = [-0.25, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../images/pdf/fig7b_G_PTCs.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
