clear all; close all; clc;

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
run_read = 'run06_initial_periodic_orbit';

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

%------------------------%
%    Add two periods     %
%------------------------%
% Periodicity
k = 2;

% Period of periodic orbit
tbp = tbp_PO / T_PO;

% Copy original periodic orbit data
xbp_plot = xbp_PO;

tbp_plot = tbp;

for i = 2: k
  xbp_plot = [xbp_plot(1:end-1, :); xbp_PO];
  tbp_plot  = [tbp_plot(1:end-1); tbp_plot(end) + tbp];
end

%%
%-------------------------------------------------------------------------%
%                          Plot: Temporal Trace                           %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Temporal Trace of Periodic Orbit';

% Figure dimensions
fig.Units = 'centimeters';
fig.Units = 'centimeters';
fig.Position = [5, 5, 8.5, 4.25];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

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

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% % Turn off all axis labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Labels     %
%---------------------%
ax.XAxis.Label.String = '$t / T_{\Gamma}$';
ax.YAxis.Label.String = '$G/Q/I$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../fig2b_periodic_orbit_temporal_trace.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
% plot2svg(filename_out);
