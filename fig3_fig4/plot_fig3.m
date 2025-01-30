clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
%----------------------------------%
%     Read Data from .mat File     %
%----------------------------------%
% Read from data structure file
load('../plot_mat_files/fig3_data.mat');

% %------------------------------------------%
% %     Read Initial Periodic Orbit Data     %
% %------------------------------------------%
% % Run string identifier
% run_PO = 'run05_initial_periodic_orbit';
% % Bifurcation data
% bd_PO = coco_bd_read(run_PO);
% 
% % Get solution label
% label_PO = coco_bd_labs(bd_PO, 'PO_PT');
% 
% % Read 'initial_PO' COLL data
% [sol_PO, data_PO] = coll_read_solution('initial_PO', run_PO, label_PO);
% 
% % State space solution
% xbp_PO = sol_PO.xbp;
% % Temporal solution
% tbp_PO = sol_PO.tbp;
% % Period
% T_PO   = sol_PO.T;
% 
% % Parameters
% p      = sol_PO.p;
% pnames = data_PO.pnames;
% 
% %-------------------------------------%
% %     Read Data: Stationary Point     %
% %-------------------------------------%
% % Read 'xpos' EP data
% [sol_pos, ~] = ep_read_solution('xpos', run_PO, label_PO);
% 
% % Stationary point
% xpos = sol_pos.x;
% 
% %--------------------------------%
% %     Read Data: Phase Reset     %
% %--------------------------------%
% % Run string indentifier
% run_PR = 'run09_phase_reset_PTC_single';
% 
% % Bifurcation data
% bd_PR = coco_bd_read(run_PR);
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
% % Get periodicity
% k = coco_bd_val(bd_PR, label_PR(1), 'k');
% 
% % Read segment 3 solution
% [sol3_run1, ~] = coll_read_solution('seg3', run_PR, label_PR(1));
% [sol3_run2, ~] = coll_read_solution('seg3', run_PR, label_PR(2));
% 
% % Read segment 4 solution
% [sol4_run1, ~] = coll_read_solution('seg4', run_PR, label_PR(1));
% [sol4_run2, ~] = coll_read_solution('seg4', run_PR, label_PR(2));
% 
% % Get segment4 state space solution
% xbp3_run1 = sol3_run1.xbp;
% xbp3_run2 = sol3_run2.xbp;
% xbp4_run1 = sol4_run1.xbp;
% xbp4_run2 = sol4_run2.xbp;
% 
% % Get segment 4 temporal solution
% tbp4_run1 = sol4_run1.tbp;
% tbp4_run2 = sol4_run2.tbp;
% 
% %----------------------------------------%
% %     Copy Unperturned Orbit k Times     %
% %----------------------------------------%
% % Time data
% tbp_normalised = tbp_PO / T_PO;
% 
% % Setup plotting data
% xbp_PO_plot = xbp_PO(:, 3);
% tbp_PO_plot = tbp_normalised;
% 
% % Copy original periodic orbit data
% for i = 2: k
%   xbp_PO_plot = [xbp_PO_plot(1:end-1); xbp_PO(:, 3)];
%   tbp_PO_plot = [tbp_PO_plot(1:end-1); tbp_PO_plot(end) + tbp_normalised];
% end
% 
% % Multiply segment time solutions
% tbp4_run1 = k * tbp4_run1;
% tbp4_run2 = k * tbp4_run2;

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
fig.Position = [5, 5, 8.5, 12.7];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(3, 2, Padding='compact', TileSpacing='compact');
ax1 = nexttile([1, 1]);
ax2 = nexttile([1, 1]);
ax3 = nexttile([1, 2]);
ax4 = nexttile([1, 2]);

% Array of axes
ax = [ax1, ax2, ax3, ax4];

% Set fontsizes
fontsize = 8;
for i = 1 : 4
  ax(i).FontSize = fontsize;
end

%----------------------------------%
%     Plot: Phase Portrait (1)     %
%----------------------------------%
% Hold axes
hold(ax1, 'on');

% Plot segment 4
plot(ax1, xbp4_run1(:, 1), xbp4_run1(:, 3), Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=1.0, DisplayName='Segment 4');

% Plot original periodic orbit
plot(ax1, xbp_PO(:, 1), xbp_PO(:, 3), Color=colours(3, :), ...
     LineWidth=2.0, DisplayName='$\Gamma$');

% Plot equilibrium point
plot(ax1, xpos(1), xpos(3), Marker='diamond', MarkerSize=5, ...
     MarkerFaceColor='r', MarkerEdgecolor='r', ...
     HandleVisibility='off');

% Marker for theta_old point
plot(ax1, xpos(1), xpos(3), Marker='o', MarkerSize=7.5, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% Hold axes
hold(ax1, 'off');

%----------------------------------%
%     Plot: Phase Portrait (2)     %
%----------------------------------%
% Hold axes
hold(ax2, 'on');

% Plot segment 4
plot(ax2, xbp4_run2(:, 1), xbp4_run2(:, 3), Color=[0.0, 0.0, 0.0, 0.5], ...
     LineWidth=1.0, DisplayName='Segment 4');

% Plot original periodic orbit
plot(ax2, xbp_PO(:, 1), xbp_PO(:, 3), Color=colours(3, :), ...
     LineWidth=2.0, DisplayName='$\Gamma$');

% Plot equilibrium point
plot(ax2, xpos(1), xpos(3), Marker='diamond', MarkerSize=5, ...
     MarkerFaceColor='r', MarkerEdgecolor='r', ...
     HandleVisibility='off');

% Marker for theta_old point
plot(ax2, xpos(1), xpos(3), Marker='o', MarkerSize=7.5, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% Hold axes
hold(ax2, 'off');

%-------------------------------%
%     Plot: Time Series (1)     %
%-------------------------------%
% Hold axes
hold(ax3, 'on');

% Plot unerperturbed orbit
max_idx = max(find(tbp_PO_plot < 12.0));
plot(ax3, tbp_PO_plot(1:max_idx), xbp_PO_plot(1:max_idx), Color=[colours(3, :)], LineWidth=1.0);

% Plot segment 4
max_idx = max(find(tbp4_run1 < 12.0));
plot(ax3, tbp4_run1(1:max_idx), xbp4_run1(1:max_idx, 3), Color=[0.0, 0.0, 0.0, 0.5], LineWidth=1.0);

% Hold axes
hold(ax3, 'off');

%-------------------------------%
%     Plot: Time Series (2)     %
%-------------------------------%
% Hold axes
hold(ax4, 'on');

% Plot unerperturbed orbit
max_idx = max(find(tbp_PO_plot < 12.0));
plot(ax4, tbp_PO_plot(1:max_idx), xbp_PO_plot(1:max_idx), Color=[colours(3, :)], LineWidth=1.0);

% Plot segment 4
max_idx = max(find(tbp4_run2 < 12.0));
plot(ax4, tbp4_run2(1:max_idx), xbp4_run2(1:max_idx, 3), Color=[0.0, 0.0, 0.0, 0.5], LineWidth=1.0);

% Hold axes
hold(ax4, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax1.XAxis.Limits = [0, 5];
ax1.YAxis.Limits = [-0.1, 20];
ax2.XAxis.Limits = ax1.XAxis.Limits;
ax2.YAxis.Limits = ax1.YAxis.Limits;

ax3.XAxis.Limits = [-0.2, 12];
ax3.YAxis.Limits = ax1.YAxis.Limits;
ax4.XAxis.Limits = ax3.XAxis.Limits;
ax4.YAxis.Limits = ax3.YAxis.Limits;

%------------------------------%
%     Axis Ticks: Settings     %
%------------------------------%
for i = 1 : 4
  ax(i).XAxis.TickDirection = 'in';
  ax(i).XAxis.MinorTick = 'on';
  ax(i).YAxis.TickDirection = 'in';
  ax(i).YAxis.MinorTick = 'on';
end

%---------------------------------%
%     Axis Ticks: ax1 and ax2     %
%---------------------------------%
% X-Axis
ax1.XAxis.TickValues = 0.0 : 1 : 5.0;
ax1.XAxis.MinorTickValues = 0.0 : 0.5 : 5.0;

% Y-Axis
ax1.YAxis.TickValues = 0.0 : 5 : 20.0;
ax1.YAxis.MinorTickValues = 0.0 : 2.5 : 20.0;

% ax2
ax2.XAxis.TickValues      = ax1.XAxis.TickValues;
ax2.XAxis.MinorTickValues = ax1.XAxis.MinorTickValues;
ax2.YAxis.TickValues      = ax1.YAxis.TickValues;
ax2.YAxis.MinorTickValues = ax1.YAxis.MinorTickValues;
ax2.XAxis.TickLabels      = ax1.XAxis.TickLabels;

%---------------------------------%
%     Axis Ticks: ax3 and ax4     %
%---------------------------------%
% X-Axis
ax3.XAxis.TickValues = 0.0 : 2.0 : 12.0;
ax3.XAxis.MinorTickValues = 0.0 : 1.0 : 12.0;

% Y-Axis
ax3.YAxis.TickValues = ax1.YAxis.TickValues;
ax3.YAxis.MinorTickValues = ax1.YAxis.MinorTickValues;

ax4.XAxis.TickValues      = ax3.XAxis.TickValues;
ax4.XAxis.MinorTickValues = ax3.XAxis.MinorTickValues;
ax4.YAxis.TickValues      = ax3.YAxis.TickValues;
ax4.YAxis.MinorTickValues = ax3.YAxis.MinorTickValues;
ax4.YAxis.TickLabels      = ax3.YAxis.TickLabels;

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% % Turn off all axis labels\
% for i = 1 : 4
%   ax(i).XAxis.TickLabels = {};
%   ax(i).YAxis.TickLabels = {};
% end

%---------------------%
%     Axis Labels     %
%---------------------%
for i = 1 : 2
  ax(i).XAxis.Label.String = '$G$';
  ax(i).YAxis.Label.String = '$I$';
end

for i = 3 : 4
  ax(i).XAxis.Label.String = '$t / T_{\Gamma}$';
  ax(i).YAxis.Label.String = '$I$';
end

%----------------------%
%     Figure Stuff     %
%----------------------%
for i = 1 : 4
  box(ax(i), 'on');
end

%---------------------%
%     Save Figure     %
%---------------------%
% filename_out = '../fig3_G_reset_phase_and_time.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
