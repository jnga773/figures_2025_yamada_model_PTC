clear all; close all; clc;

%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
%----------------------------------%
%     Read Data from .mat File     %
%----------------------------------%
load('../plot_mat_files/fig4_data.mat');

%----------------------------------%
%     Read Data from COCO Data     %
%----------------------------------%
% % Bifurcation data
% bd_PR = coco_bd_read(run_in);
% 
% % Get solution labels
% label_PR = coco_bd_labs(bd_PR, 'SP');
% 
% % Get theta_old values
% theta_old_run1 = coco_bd_val(bd_PR, label_PR(1), 'theta_old');
% theta_old_run2 = coco_bd_val(bd_PR, label_PR(2), 'theta_old');
% 
% % Get A_perturb value
% A_perturb = coco_bd_val(bd_PR, label_PR(1), 'A_perturb');
% 
% % Read theta_old and theta_new values from bifurcation data
% theta_old_plot = coco_bd_col(bd_PR, 'theta_old');
% theta_new_plot = coco_bd_col(bd_PR, 'theta_new');

%%
%-------------------------------------------------------------------------%
%                       Plot Phase Transition Curve                       %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(2); clf;
fig.Name = 'Single PTC';

% Figure dimensions
fig.Units = 'centimeters';
fig.Position = [5, 5, 6.0, 6.0];

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
% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='-', Color=colours(3, :), LineWidth=1.5, ...
     HandleVisibility='off');

% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

% Grey lines at theta_old = 0 and 0.3
xline(ax, theta_old_run1, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);
xline(ax, theta_old_run2, LineStyle='-', Color=[0, 0, 0, 0.5], LineWidth=1);

%-------------------%
%     Plot: PTC     %
%-------------------%
% Plot PTC: run1
plot(ax, theta_old_plot, theta_new_plot, Color='k', LineStyle='-', LineWidth=1.5);
plot(ax, theta_old_plot, theta_new_plot+1, Color='k', LineStyle='-', LineWidth=1.5);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-0.0, 1.0];
ax.YAxis.Limits = [-0.1, 1.0];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = -1.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -1.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% % Turn off all axis labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$\theta_{\mathrm{o}}$';
ax.YAxis.Label.String = '$\theta_{\mathrm{n}}$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
% filename_out = '../fig4_G_PTC_single.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
