clear all; close all; clc;

%%
%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Read from data structure file
load('./data_files/fig3_data.mat');
load('./data_files/initial_PO.mat');

%-------------------------%
%     Read Parameters     %
%-------------------------%
% Periodicity
% k = p_run1(6);
% Period
% T_PO = p_run1(5);

% Read A_perturb
% A_perturb_run1 = p_run1(11);
% A_perturb_run2 = p_run2(11);

fprintf('A_perturb(1) = %.4f\n', A_perturb);
fprintf('A_perturb(2) = %.4f\n\n', A_perturb);

% Read theta_old
% theta_old_run1 = p_run1(7);
% theta_old_run2 = p_run2(7);

fprintf('theta_old(1) = %.4f\n', theta_old_run1);
fprintf('theta_old(2) = %.4f\n\n', theta_old_run2);

%----------------------------------------%
%     Copy Unperturned Orbit k Times     %
%----------------------------------------%
% Time data
tbp = tbp_PO / T_PO;

% Copy original periodic orbit data
X_PO_plot = xbp_PO(:, 3);
tbp_plot  = tbp;
for i = 2: k
  X_PO_plot = [X_PO_plot(1:end-1); xbp_PO(:, 3)];
  tbp_plot  = [tbp_plot(1:end-1); tbp_plot(end) + tbp];
end

% Multiply segment time solutions
tbp4_run1 = k * tbp4_run1;
tbp4_run2 = k * tbp4_run2;

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Default line colours
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
fig.Units = 'inches'; fig.Position = [3, 3, 3.3, 5]; fig.PaperSize = [3.3, 5];

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
plot(ax1, xbp3_run1(1, 1), xbp3_run1(1, 3), Marker='o', MarkerFaceColor='k', ...
     MarkerEdgeColor='k');

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
plot(ax2, xbp3_run2(1, 1), xbp3_run2(1, 3), Marker='o', MarkerFaceColor='k', ...
     MarkerEdgeColor='k');

% Hold axes
hold(ax2, 'off');

%-------------------------------%
%     Plot: Time Series (1)     %
%-------------------------------%
% Hold axes
hold(ax3, 'on');

% Plot unerperturbed orbit
max_idx = max(find(tbp_plot < 12.0));
plot(ax3, tbp_plot(1:max_idx), X_PO_plot(1:max_idx), Color=[colours(3, :)], LineWidth=1.0);

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
max_idx = max(find(tbp_plot < 12.0));
plot(ax4, tbp_plot(1:max_idx), X_PO_plot(1:max_idx), Color=[colours(3, :)], LineWidth=1.0);

% Plot segment 4
max_idx = max(find(tbp4_run2 < 12.0));
plot(ax4, tbp4_run2(1:max_idx), xbp4_run2(1:max_idx, 3), Color=[0.0, 0.0, 0.0, 0.5], LineWidth=1.0);

% Hold axes
hold(ax4, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax1.XAxis.Limits = [0, 5];
ax1.YAxis.Limits = [-1, 20];
ax2.XAxis.Limits = ax1.XAxis.Limits;
ax2.YAxis.Limits = ax1.YAxis.Limits;

ax3.XAxis.Limits = [-0.2, 12.2];
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

% Labels
ax2.YAxis.TickLabels = {};

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

% Labels
ax3.XAxis.TickLabels = {};

ax4.XAxis.TickValues      = ax3.XAxis.TickValues;
ax4.XAxis.MinorTickValues = ax3.XAxis.MinorTickValues;
ax4.YAxis.TickValues      = ax3.YAxis.TickValues;
ax4.YAxis.MinorTickValues = ax3.YAxis.MinorTickValues;
ax4.YAxis.TickLabels      = ax3.YAxis.TickLabels;

%--------------------------%
%     Axis Tick Labels     %
%--------------------------%
% Turn off all axis labels\
for i = 1 : 4
  ax(i).XAxis.TickLabels = {};
  ax(i).YAxis.TickLabels = {};
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
filename_out = '../images/pdf/fig3_G_reset_phase_and_time.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
% plot2svg('../images/fig2/testsvg.svg');
% fig2svg('../images/fig2/testsvg.svg', fig, [3.3, 5]);



