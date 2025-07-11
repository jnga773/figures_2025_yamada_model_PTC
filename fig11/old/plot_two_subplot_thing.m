close all; clear;

%%
%-------------------------------------------------------------------------%
%                              Calculate Data                             %
%-------------------------------------------------------------------------%
% Calculate data
[theta_old, theta_perturb, A_perturb] = calc_intersection_points();

% Normalise theta_perturb by pi
theta_perturb = theta_perturb / pi;

%-----------------------------------------%
%     Sort Data: Increasing theta_old     %
%-----------------------------------------%
% Sort
[~, sort_idx] = sort(theta_old);

% Sort by theta_perturb
theta_old     = theta_old(sort_idx);
theta_perturb = theta_perturb(sort_idx);
A_perturb     = A_perturb(sort_idx);

%---------------------------------------------%
%     Sort Data: Increasing theta_perturb     %
%---------------------------------------------%
% % Sort
% [~, sort_idx] = sort(theta_perturb);
% 
% % Sort by theta_perturb
% theta_old     = theta_old(sort_idx);
% theta_perturb = theta_perturb(sort_idx);
% A_perturb     = A_perturb(sort_idx);

%------------------------------------------%
%     Shift Data: By max theta_perturb     %
%------------------------------------------%
% Find max point of theta_old
[max_val, max_idx] = max(theta_perturb);

% Shift data around
theta_perturb = [theta_perturb(1:max_idx)-2; theta_perturb(max_idx+1:end)];
theta_old     = [theta_old(1:max_idx); theta_old(max_idx+1:end)];
A_perturb     = [A_perturb(1:max_idx); A_perturb(max_idx+1:end)];

%---------------------%
%     Extend Data     %
%---------------------%
theta_perturb = [theta_perturb-2; theta_perturb; theta_perturb+2];
theta_old     = [theta_old-1; theta_old; theta_old+1];
A_perturb     = [A_perturb; A_perturb; A_perturb];

%%
%--------------------------------------------%
%     Interpolate Data: By theta_perturb     %
%--------------------------------------------%
% Ranges
range_min = -1.0;
range_max = 2.5;
drange    = 1e-3;

% Get data over range theta_perturb in [-0.5, 1.0]
range_idx = theta_perturb >= range_min & theta_perturb <= range_max;

% Make data smaller range
theta_perturb_small_range = theta_perturb(range_idx);
theta_old_small_range     = theta_old(range_idx);
A_perturb_small_range     = A_perturb(range_idx);

[~, unique_idx] = unique(theta_perturb_small_range);
theta_perturb_small_range = theta_perturb_small_range(unique_idx);
theta_old_small_range = theta_old_small_range(unique_idx);
A_perturb_small_range = A_perturb_small_range(unique_idx);

% Smooth theta_perturb data to interpolate over
theta_perturb_interpolate = range_min : drange : range_max;
theta_perturb_interpolate = theta_perturb_interpolate';
% theta_perturb_interpolate = unique(theta_perturb_interpolate);

% Inteprolate data
theta_old_interpolate = interp1(theta_perturb_small_range, theta_old_small_range, theta_perturb_interpolate);
A_perturb_interpolate = interp1(theta_perturb_small_range, A_perturb_small_range, theta_perturb_interpolate);

% Halve theta_perturb
theta_perturb_interpolate = 0.5 * theta_perturb_interpolate;

% Print values at theta_perturb = 0 and "pi"
idx_G    = find(theta_perturb_interpolate == 0.0);
idx_I    = find(theta_perturb_interpolate == 0.25);
% idx_half = find(theta_perturb_interpolate == 0.125);

perturb_angle = 45;
idx_half = find(round(theta_perturb_interpolate, 3) == round(deg2rad(perturb_angle) / (2 * pi), 3));
idx_half = idx_half(1);

% Intersection for G perturbation
fprintf('G intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
        mod(theta_old_interpolate(idx_G), 1), A_perturb_interpolate(idx_G));

% Intersection for I perturbation
fprintf('I intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
        mod(theta_old_interpolate(idx_I), 1), A_perturb_interpolate(idx_I));

% Intersection for 45deg perturbation
fprintf('45deg intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
        mod(theta_old_interpolate(idx_half), 1), A_perturb_interpolate(idx_half));

%%
%-------------------------------------------------------------------------%
%                         Plot Data: 2D Comparison                        %
%-------------------------------------------------------------------------%
colours = colororder();

% Setup figure
fig = figure(2); clf;
fig.Name = 'Theta_old vs theta_perturb vs A_perturb';

% Figure dimensions
fig.Units = 'inches';
fig.Position = [4, 4, 8, 8];

% Axis setup: Manual padding
% ax = gca();
% ax.Position = [0.01, 0.01, 0.98, 0.98];

% Axis setup: Tiled layout
tiles = tiledlayout(2, 1, Padding='compact', TileSpacing='tight');
ax1 = nexttile;
ax2 = nexttile;
ax = [ax1, ax2];

% Fontsize
for i = 1 : 2
  ax(i).FontSize = 12;
end

set(ax, 'TickLabelInterpreter','latex');  
set(ax, 'DefaultTextInterpreter', 'latex');

%------------------------------------------%
%     Plot: theta_perturb vs theta_old     %
%------------------------------------------%
hold(ax(1), 'on');

% Highlight area over 0.0 <= theta_perturb <= 0.5 pi
patch(ax(1), [0, 0, 0.25, 0.25], [0, 1, 1, 0], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

% Plot interpolated data
plot(ax(1), theta_perturb_interpolate, theta_old_interpolate, ...
     Color=colours(1, :), LineStyle='-', LineWidth=2.5);
plot(ax(1), theta_perturb_interpolate+2, theta_old_interpolate, ...
     Color=colours(1, :), LineStyle='-', LineWidth=2.5);
plot(ax(1), theta_perturb_interpolate, theta_old_interpolate-1, ...
     Color=colours(1, :), LineStyle='-', LineWidth=2.5);

% Add intersection points
plot(ax(1), theta_perturb_interpolate(idx_G), theta_old_interpolate(idx_G), ...
     LineStyle='none', Marker='o', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');
plot(ax(1), theta_perturb_interpolate(idx_I), theta_old_interpolate(idx_I), ...
     LineStyle='none', Marker='diamond', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');
plot(ax(1), theta_perturb_interpolate(idx_half), theta_old_interpolate(idx_half), ...
     LineStyle='none', Marker='square', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');

% Add text
plot_text = sprintf('$\\vartheta_{\\mathrm{o}}^{(0^{\\circ})} = %.6f$', theta_old_interpolate(idx_G));
text(ax(1), theta_perturb_interpolate(idx_G), theta_old_interpolate(idx_G)-0.05, ...
     plot_text);
plot_text = sprintf('$\\vartheta_{\\mathrm{o}}^{(90^{\\circ})} = %.6f$', theta_old_interpolate(idx_I));
text(ax(1), theta_perturb_interpolate(idx_I)+0.015, theta_old_interpolate(idx_I), ...
     plot_text);
plot_text = sprintf('$\\vartheta_{\\mathrm{o}}^{(45^{\\circ})} = %.6f$', theta_old_interpolate(idx_half));
text(ax(1), theta_perturb_interpolate(idx_half)-0.1, theta_old_interpolate(idx_half)+0.05, ...
     plot_text);

hold(ax(1), 'off')

%------------------------------------------%
%     Plot: theta_perturb vs A_perturb     %
%------------------------------------------%
hold(ax(2), 'on');

% Highlight area over 0.0 <= theta_perturb <= 0.5 pi
patch(ax(2), [0, 0, 0.25, 0.25], [0, 20, 20, 0], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

% Plot interpolated data
plot(ax(2), theta_perturb_interpolate, A_perturb_interpolate, ...
     Color=colours(1, :), LineStyle='-', LineWidth=2.5);
plot(ax(2), theta_perturb_interpolate+2, A_perturb_interpolate, ...
     Color=colours(1, :), LineStyle='-', LineWidth=2.5);

% Add intersection points
plot(ax(2), theta_perturb_interpolate(idx_G), A_perturb_interpolate(idx_G), ...
     LineStyle='none', Marker='o', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');
plot(ax(2), theta_perturb_interpolate(idx_I), A_perturb_interpolate(idx_I), ...
     LineStyle='none', Marker='diamond', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');
plot(ax(2), theta_perturb_interpolate(idx_half), A_perturb_interpolate(idx_half), ...
     LineStyle='none', Marker='square', MarkerSize=5, ...
     MarkerFaceColor='k', MarkerEdgeColor='k');

% Add text
plot_text = sprintf('$A^{(0^{\\circ})} = %.6f$', A_perturb_interpolate(idx_G));
text(ax(2), theta_perturb_interpolate(idx_G)-0.1, A_perturb_interpolate(idx_G)+0.75, ...
     plot_text);
plot_text = sprintf('$A^{(90^{\\circ})} = %.6f$', A_perturb_interpolate(idx_I));
text(ax(2), theta_perturb_interpolate(idx_I)-0.175, A_perturb_interpolate(idx_I), ...
     plot_text);
plot_text = sprintf('$A^{(45^{\\circ})} = %.6f$', A_perturb_interpolate(idx_half));
text(ax(2), theta_perturb_interpolate(idx_half)-0.1, A_perturb_interpolate(idx_half)+1.5, ...
     plot_text);

hold(ax(2), 'off')

%----------------------------%
%     Axis Ticks: Axis 1     %
%----------------------------%
% X-Axis: theta_perturb
ax(1).XAxis.MinorTick = 'on';
ax(1).XAxis.TickValues = -0.25 : 0.25 : 1.5;
ax(1).XAxis.MinorTickValues = -0.25 : 0.125 : 1.5;
ax(1).XAxis.TickLabels = {};

% Y-Axis: theta_old
ax(1).YAxis.MinorTick = 'on';
ax(1).YAxis.TickValues = 0.0 : 0.2 : 1.0;
ax(1).YAxis.MinorTickValues = 0.0 : 0.05 : 1.0;

%----------------------------%
%     Axis Ticks: Axis 2     %
%----------------------------%
% X-Axis: theta_perturb
ax(2).XAxis.MinorTick = ax(1).XAxis.MinorTick;
ax(2).XAxis.TickValues = ax(1).XAxis.TickValues;
ax(2).XAxis.MinorTickValues = ax(1).XAxis.MinorTickValues;
% ax(2).XAxis.TickLabels = {'$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', ...
%                           '$\pi$', '$\frac{3 \pi}{2}$', '$2 \pi$', '$\frac{5 \pi}{2}$'};
% ax(2).XAxis.TickLabels = {'$-\pi/2$', '$0$', '$\pi/2$', ...
%                           '$\pi$', '$3\pi/2$', '$2 \pi$', '$5\pi/2$'};


% Y-Axis: A_perturb
ax(2).YAxis.TickDirection = 'in';
ax(2).YAxis.MinorTick = 'on';
ax(2).YAxis.TickValues = 0.0 : 3.0 : 15.0;
ax(2).YAxis.MinorTickValues = 0.0 : 1.0 : 15.0;

%---------------------%
%     Axis Limits     %
%---------------------%
% X-Axis: theta_perturb
% ax(1).XAxis.Limits = [-0.5, 2.5];
ax(1).XAxis.Limits = [-0.125, 1.0];
ax(2).XAxis.Limits = ax(1).XAxis.Limits;

% Y-Axis: theta_old
ax(1).YAxis.Limits = [0.0, 1.0];

% Y-Axis: A_perturb
ax(2).YAxis.Limits = [0.0, 15.0];

%---------------------%
%     Axis Labels     %
%---------------------%
% X-Axis: theta_perturb
ax(2).XAxis.Label.String = '$\varphi_{\mathrm{d}} / 2 \pi$';

% Y-Axis: theta_old
ax(1).YAxis.Label.String = '$\vartheta_{\mathrm{o}}$';

% Y-Axis: A_perturb
ax(2).YAxis.Label.String = '$A$';

%----------------------%
%     Figure Stuff     %
%----------------------%
for i = 1 : 2
  box(ax(i), 'on');
  % grid(ax(i), 'on');
end

% Save figure
filename_out = './theta_old_A_perturb_intersection.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
