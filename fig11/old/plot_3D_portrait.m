close all; clear;

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

%%
%-------------------------------------------------------------------------%
%                          Plot Data: 3D Portrait                         %
%-------------------------------------------------------------------------%
colours = colororder();

% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';

% Figure dimensions
fig.Units = 'inches';
fig.Position = [3, 3, 12, 8];

% Axis setup: Manual padding
% ax = gca();
% ax.Position = [0.01, 0.01, 0.98, 0.98];

% Axis setup: Tiled layout
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

% Fontsize
ax.FontSize = 9;

%--------------%
%     Plot     %
%--------------%
hold(ax, 'on');

% Plot original periodic orbit
plot3(ax, xbp_PO(:, 1), xbp_PO(:, 2), xbp_PO(:, 3), Color=colours(3, :), ...
      LineWidth=2.0);

% Plot stable manifold
plot3(ax, Wq_s(:, 1), Wq_s(:, 2), Wq_s(:, 3), ...
      Color=colours(1, :), LineWidth=2.0);

% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), Marker='o', MarkerSize=7.5, LineStyle='none', ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

% % Cycle through all points and plot perturbation vector towards stable
% % manifold.
% for i = 1 : length(Q_plot(:, 1))
%   if theta_perturb(i) <= 0.5*pi
%     % Plot
%     plot3(ax, G_plot(i, :), Q_plot(i, :), I_plot(i, :), Color=colours(4, :), ...
%           HandleVisibility='off');
%   end
% end

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Get indices for theta_perturb <= 0.5 pi
theta_idx = theta_perturb <= 0.5 * pi;

% Get surface data
G_surf = G_plot(theta_idx, :);
Q_surf = Q_plot(theta_idx, :);
I_surf = I_plot(theta_idx, :);

% Plot surface
surf(ax, G_surf, Q_surf, I_surf, ...
     EdgeColor='none', FaceColor=colours(4, :), FaceAlpha=0.5);

% Plot specific ones
[~, min_idx] = min(theta_perturb(theta_idx));
[~, max_idx] = max(theta_perturb(theta_idx));
[~, mid_idx] = min(abs(theta_perturb(theta_idx) - 0.25 * pi));
theta_idx_2 = [min_idx, mid_idx, max_idx];

for i = 1 : length(theta_idx_2)
  idx = theta_idx_2(i);
  plot3(ax, G_surf(idx, :), Q_surf(idx, :), I_surf(idx, :), ...
        LineStyle='-', Color=colours(4, :), LineWidth=4.0);
end

hold(ax, 'off')
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
view(45, 6.0);