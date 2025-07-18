close all; clear; clc;

%%
%-------------------------------------------------------------------------%
%                              Calculate Data                             %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s', 'xpos');

% Calculate data
[theta_old, theta_perturb, A_perturb] = calc_intersection_points();

% Get perturbation vector
d_vec = A_perturb .* [cos(theta_perturb), sin(theta_perturb)];

% Plotting vector
G_plot = [xbp_PO(:, 1), xbp_PO(:, 1) + d_vec(:, 1)];
Q_plot = [xbp_PO(:, 2), xbp_PO(:, 2)];
I_plot = [xbp_PO(:, 3), xbp_PO(:, 3) + d_vec(:, 2)];

%----------------------%
%     Plot Colours     %
%----------------------%
% Default colour order
colours = colororder();

%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
width = 3.0;
height = 1.6;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%--------------%
%     Plot     %
%--------------%
% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), ...
      Color=colours(3, :), ...
      LineWidth=2.0);

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), ...
      LineWidth=2.0);

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Linewidth
lw = 3.5;
colour = colours(7, :);

% Get indices for theta_perturb <= 0.5 pi
theta_idx = theta_perturb <= 0.5 * pi;

% Get surface data
G_surf = G_plot(theta_idx, :);
Q_surf = Q_plot(theta_idx, :);
I_surf = I_plot(theta_idx, :);

% Plot surface
surf(ax, G_surf, Q_surf, I_surf, ...
     EdgeColor='none', FaceColor=colour, FaceAlpha=0.25);

% Plot specific ones
[~, min_idx] = min(theta_perturb(theta_idx));
[~, max_idx] = max(theta_perturb(theta_idx));
[~, mid_idx] = min(abs(theta_perturb(theta_idx) - 0.25 * pi));
% theta_idx_2 = [min_idx, mid_idx, max_idx];
theta_idx_2 = [min_idx, max_idx];

for i = 1 : length(theta_idx_2)
  idx = theta_idx_2(i);
  plot3(ax, G_surf(idx, :), Q_surf(idx, :), I_surf(idx, :), ...
        LineStyle='-', Color=colour, LineWidth=lw);
end

plot3(ax, G_surf(:, 1), Q_surf(:, 1), I_surf(:, 1), ...
      LineStyle='-', Color=colour, LineWidth=lw);
plot3(ax, G_surf(:, 2), Q_surf(:, 2), I_surf(:, 2), ...
      LineStyle='-', Color=colour, LineWidth=lw);

%---------------------------------%
%     Plot: Equilibrium Point     %
%---------------------------------%
% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), ...
      Marker='o', MarkerSize=4, ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 3.0];
ax.YAxis.Limits = [0.0, 3.0];
ax.ZAxis.Limits = [0.0, 5];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = [];

% Y-Axis
ax.YAxis.TickValues = [];

% Z-Axis
ax.ZAxis.TickValues = [];

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$G$');
% ylabel(ax, '$Q$');
% zlabel(ax, '$I$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};
ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

% Grid lines
% ax.GridLineWidth = 0.5;
% ax.GridColor = 'black';
% ax.GridAlpha = 0.25;

% 3D plot view
view(45, 6.0);
% view(-45, 6);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig11a_inset1.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
