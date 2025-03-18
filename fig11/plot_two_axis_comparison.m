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
theta_old     = [theta_old(1:max_idx); theta_old(max_idx+1:end)];
theta_perturb = [theta_perturb(1:max_idx)-2; theta_perturb(max_idx+1:end)];
A_perturb     = [A_perturb(1:max_idx); A_perturb(max_idx+1:end)];

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
fig.Position = [4, 4, 8, 5];

% Axis setup: Manual padding
% ax = gca();
% ax.Position = [0.01, 0.01, 0.98, 0.98];

% Axis setup: Tiled layout
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;

% Fontsize
ax.FontSize = 12;

set(ax, 'TickLabelInterpreter','latex');  
set(ax, 'DefaultTextInterpreter', 'latex');

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%--------------------------------------------%
%     Plot: Highlight theta_perturb Area     %
%--------------------------------------------%
% Highlight area over 0.0 <= theta_perturb <= 0.5 pi
patch(ax(1), [0, 0, 0.5, 0.5], [0, 1, 1, 0], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

%------------------------------------------%
%     Plot: theta_perturb vs theta_old     %
%------------------------------------------%
% Left axis
yyaxis(ax, 'left');

% Plot: theta_old vs theta_perturb
plot(ax(1), theta_perturb, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
plot(ax(1), theta_perturb+2, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);

%------------------------------------------%
%     Plot: theta_perturb vs A_perturb     %
%------------------------------------------%
% Right axis
yyaxis(ax, 'right');

% Plot: theta_old vs theta_perturb
plot(ax, theta_perturb, A_perturb, Color=colours(2, :), LineStyle='-', LineWidth=2.5);
plot(ax, theta_perturb+2, A_perturb, Color=colours(2, :), LineStyle='-', LineWidth=2.5);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis: theta_perturb
ax.XAxis.MinorTick = 'on';
ax.XAxis.TickValues = -0.5 : 0.5 : 2.5;
ax.XAxis.MinorTickValues = -0.5 : 0.25 : 2.5;
% ax.XAxis.TickLabels = {'$-\frac{\pi}{2}$', '$0$', '$\frac{\pi}{2}$', ...
%                        '$\pi$', '$\frac{3 \pi}{2}$', '$2 \pi$', '$\frac{5 \pi}{2}$'};
ax.XAxis.TickLabels = {'$-\pi/2$', '$0$', '$\pi/2$', ...
                       '$\pi$', '$3\pi/2$', '$2 \pi$', '$5\pi/2$'};

% Y-Axis 1: theta_old
ax.YAxis(1).MinorTick = 'on';
ax.YAxis(1).TickValues = 0.0 : 0.1 : 1.0;
ax.YAxis(1).MinorTickValues = 0.0 : 0.05 : 1.0;

% Y-Axis 2: A_perturb
ax.YAxis(2).MinorTick = 'on';
ax.YAxis(2).TickValues = 0.0 : 3.0 : 15.0;
ax.YAxis(2).MinorTickValues = 0.0 : 1.0 : 15.0;

%---------------------%
%     Axis Limits     %
%---------------------%
% X-Axis: theta_perturb
ax.XAxis.Limits = [-0.5, 2.5];

% Y-Axis 1: theta_old
ax.YAxis(1).Limits = [0.0, 1.0];

% Y-Axis 2: A_perturb
ax.YAxis(2).Limits = [0.0, 15.0];

%---------------------%
%     Axis Labels     %
%---------------------%
% X-Axis: theta_perturb
ax.XAxis.Label.String = '$\varphi_{\mathrm{p}}$';

% Y-Axis 1: theta_perturb
ax.YAxis(1).Label.String = '$\vartheta_{\mathrm{o}}$';

% Y-Axis 2: A_perturb
ax.YAxis(2).Label.String = '$A$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');
