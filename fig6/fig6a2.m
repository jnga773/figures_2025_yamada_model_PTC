clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read data from initial periodic orbit data .mat file
filename_data = './data_files/initial_PO.mat';
% Load data
load(filename_data);

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';

% Figure dimensions
fig.Units = 'centimeters';
fig.Position = [3, 3, 8.5, 6];
% fig.Position = [3, 3, 7.083, 5];
% fig.Position = [3, 3, 6.375, 4.5];
% fig.Position = [3, 3, 5.667, 4];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='none');
ax = nexttile;
ax.FontSize = 8;

% set(fig, 'DefaultFigureRenderer', 'painters');
% set(fig,'renderer','Painters');

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%--------------%
%     Plot     %
%--------------%
% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), Color=colours(3, :), ...
      LineWidth=2.0);

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), Marker='o', MarkerSize=7.5, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%----------------------------%
%     Plot: Perturbation     $
%----------------------------%
% List of perturbations
A_p = [0.55, 1.0, 1.4, 1.8];

plot_colours = {colours(6, :), colours(7, :), colours(8, :), colours(9, :)};
lw = 1.0;

for i = 1 : length(A_p)
  plot3(ax, smooth(xbp_PO(:, 1)+A_p(i)), smooth(xbp_PO(:, 2)), smooth(xbp_PO(:, 3)), ...
        Color=plot_colours{i}, LineStyle='-', LineWidth=lw)%, ...
        % Marker='none', MarkerFaceColor=colour(1:3), MarkerEdgeColor=colour(1:3));
end

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 6.0];
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, 20];


%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 2.0 : 6.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 6.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 2.0 : 4.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 1.0 : 4.0;

% Z-Axis
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickValues = 0.0 : 5.0 : 20.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 0.0 : 2.5 : 20.0;

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% Turn off all axis labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};
ax.ZAxis.TickLabels = {};

%---------------------%
%     Axis Labels     %
%---------------------%
% ax.XAxis.Label.String = '$G$';
% ax.YAxis.Label.String = '$Q$';
% ax.ZAxis.Label.String = '$I$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% 3D plot view
% view(45, 10.0);
view(45, 6);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../images/pdf/fig7a_phase_portrait.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
% saveas(fig, filename_out, 'svg');
% exportgraphics(fig, filename_out, 'ContentType', 'vector', 'Units', 'centimeters')
% export_fig(filename_out, '-pdf', fig);
% plot2svg(filename_out);
