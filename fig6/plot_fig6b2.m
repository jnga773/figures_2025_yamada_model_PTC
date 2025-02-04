close all; clear all; clc

%-------------------------------------------------------------------------%
%%                               Read Data                               %%
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig6_data.mat');

% Data indices to plot
% plot_idx = 1:4;
plot_idx = 4:7;

% Colours
% Green     (#2ca02c) = [ 44, 160,  44] ./ 255
% Yellowg   (#7fbf08) = [127, 191,   8] ./ 255
% Chartreus (#bcbd22) = [188, 189,  34] ./ 255
% Orange    (#ff7f0e) = [255, 126,  14] ./ 255
% Red       (#d62728) = [214,  39,  40] ./ 255
% Pink      (#e38ab7) = [227, 138, 183] ./ 255
% Purple    (#9467bd) = [148, 103, 189] ./ 255
% Blue      (#4649a6) = [ 70,  73, 166] ./ 255

% Plot colours
plot_colours = {[127, 191,   8] ./ 255;
                [188, 189,  34] ./ 255;
                [255, 126,  14] ./ 255;
                [214,  39,  40] ./ 255;
                [227, 138, 183] ./ 255;
                [148, 103, 189] ./ 255;
                [ 70,  73, 166] ./ 255};

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(1); clf;
fig.Name = 'PTC Scan';

% Figure dimensions
% fig.Units = 'centimeters';
fig.Units = 'inches';
fig.Position = [5, 5, 4.5, 9];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% % Axis setup: Manual padding
% ax = gca();
% ax.Position = [0.01, 0.01, 0.98, 0.98];

% Axis setup: Tiled layout
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

% Fontsize
ax.FontSize = 9;

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

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
ax.DataAspectRatio = [1, 1, 1];

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = -0.5 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -0.5 : 0.25 : 1.0;

ax.YAxis.TickValues = -0.5 : 0.5 : 2.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.5 : 0.25 : 2.5;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
xlabel(ax, '$\theta_{\mathrm{o}}$');
ylabel(ax, '$\theta_{\mathrm{n}}$');

% % Turn off all tick labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1.0];
ax.YAxis.Limits = [-0.25, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
% filename_out = '../fig6b2_G_PTCs.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
