clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read from .mat file
% % Read data from initial periodic orbit data .mat file
% filename_data = '../plot_mat_files/fig2_data.mat';
% % Load data
% load(filename_data);

% Read from COCO data
% Run string identifier
run_read = 'run09_stable_manifold_close_eps';

% Read bd file
bd_read = coco_bd_read(run_read);

% Solution label to plot (take one less just in case MX)
label_read = max(coco_bd_labs(bd_read, '')) - 1;

%-----------------------------------%
%     Read Data: Periodic Orbit     %
%-----------------------------------%
% Read COCO solution
[sol_PO, data_PO] = coll_read_solution('PO_stable', run_read, label_read);

% State space solution
xbp_PO = sol_PO.xbp;
% Temporal solution
tbp_PO = sol_PO.tbp;
% Period
T_PO   = sol_PO.T;

% Parameters
p      = sol_PO.p;
pnames = data_PO.pnames;

%--------------------------------------%
%     Read Data: Stationary Points     %
%--------------------------------------%
% Read COCO solutions
[sol_0, ~]   = ep_read_solution('x0', run_read, label_read);
[sol_pos, ~] = ep_read_solution('xpos', run_read, label_read);
[sol_neg, ~] = ep_read_solution('xneg', run_read, label_read);

% State space solutions
x0   = sol_0.x;
xpos = sol_pos.x;
xneg = sol_neg.x;

%-----------------------------%
%     Read Data: Manifold     %
%-----------------------------%
% Read stable manifold solutions
[sol1, ~] = coll_read_solution('W1', run_read, label_read);
[sol2, ~] = coll_read_solution('W2', run_read, label_read);

% State space solutions
W1 = sol1.xbp;
W2 = sol2.xbp;

% Append to single array
Wq_s = [W1; flip(W2)];

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
fig.Units = 'centimeters'; fig.Position = [5, 5, 6, 6]; fig.PaperSize = fig.Position(3:4);

tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

%--------------%
%     Plot     %
%--------------%
% Hold axes
hold(ax, 'on');

% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), Color=colours(3, :), ...
      LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), Marker='diamond', MarkerSize=10, ...
      MarkerFaceColor='r', MarkerEdgecolor='r');

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), LineWidth=2.0);

% Hold axes
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.YAxis.Limits = [0.0, 4.0];
ax.ZAxis.Limits = [0.0, 21];

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
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};
% ax.ZAxis.TickLabels = {};

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$G$';
ax.YAxis.Label.String = '$Q$';
ax.ZAxis.Label.String = '$I$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5; ax.GridColor = 'black'; ax.GridAlpha = 0.25;

% 3D plot view
view(45, 10.0);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../fig2a_phase_portrait.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
