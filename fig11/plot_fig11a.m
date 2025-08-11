%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig11_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Default colour order
colours = colororder();

%-------------------%
%     Sort Data     %
%-------------------%
% Normalise tbp_PO
tbp = tbp_PO / T_PO;

% Indiices for theta_perturb between 0 and 0.25
idx_P = (theta_perturb >= 0.0) & (theta_perturb <= 0.25);
% Find new theta_old
theta_old_sort = theta_old(idx_P);

% Get indices for periodic orbit data that aren't in theta_gamma
idxs1 = tbp <= min(theta_old_sort);
idxs2 = tbp >= max(theta_old_sort);

% Get two plotting points
xbp_PO1 = xbp_PO(idxs1, :);
xbp_PO2 = xbp_PO(idxs2, :);

% Indices for Wsq outside of xbp_Wsq
idxs1 = Wq_s(:, 2) <= min(xbp_Wsq(:, 2));
idxs2 = Wq_s(:, 2) >= max(xbp_Wsq(:, 2));

% Get two plotting points
Wqs1 = Wq_s(idxs1, :);
Wqs2 = Wq_s(idxs2, :);

% Plot circles for DTC
A_DTC = [0.1, 0.724236, 25];
theta_old_DTC = 0.339386;

% Find tbp index matching theta_old_DTC
idx_DTC = find(round(tbp, 3) == round(theta_old_DTC, 3));

gamma_DTC = xbp_PO(idx_DTC, :);

% Create circles
theta_plot = 0 : 0.01 : 2*pi;
theta_plot = theta_plot';

DTCs = {};
for idx = 1 : length(A_DTC)
  DTCs{idx} = gamma_DTC + A_DTC(idx) * [cos(theta_plot), 0 * theta_plot, sin(theta_plot)];
end


%%
%-------------------------------------------------------------------------%
%                         Plot: 3D Phase Portrait                         %
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
% width = 8.0;
% height = 5.0;
width = 15.0;
height = 5.0;

% Set figure size
set_figure_dimensions(width, height, scale=1);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%------------------------------%
%     Plot: Periodic Orbit     %
%------------------------------%
% Linewidth
lw = 2.0;
% Special colour
special_colour = colours(7, :);
% Alpha transparency
alpha = 1.0;

% Plot original periodic orbit
plot3(ax, xbp_PO1(:, 1), xbp_PO1(:, 2), xbp_PO1(:, 3), ...
      Color=[colours(3, 1:3), alpha], ...
      LineWidth=lw);
plot3(ax, xbp_PO2(:, 1), xbp_PO2(:, 2), xbp_PO2(:, 3), ...
      Color=[colours(3, 1:3), alpha], ...
      LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_gamma(:, 1), xbp_gamma(:, 2), xbp_gamma(:, 3), ...
      Color='g', LineWidth=lw, LineStyle='-');

%-------------------------------%
%     Plot: Stable Manifold     %
%-------------------------------%
% Plot original stable manifold
plot3(ax, Wqs1(:, 1), Wqs1(:, 2), Wqs1(:, 3), ...
      Color=[colours(1, 1:3), alpha], ...
      LineWidth=lw);
plot3(ax, Wqs2(:, 1), Wqs2(:, 2), Wqs2(:, 3), ...
      Color=[colours(1, 1:3), alpha], ...
      LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 2), xbp_Wsq(:, 3), ...
      Color=colours(10, :), LineWidth=lw, LineStyle='-');

%---------------------------------%
%     Plot: Equilibrium Point     %
%---------------------------------%
% Plot equilibrium point
plot3(ax, xpos(1), xpos(2), xpos(3), ...
      Marker='o', MarkerSize=4, ...
      MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Create surface plots
G_surf = [xbp_gamma(:, 1), xbp_Wsq(:, 1)];
Q_surf = [xbp_gamma(:, 2), xbp_Wsq(:, 2)];
I_surf = [xbp_gamma(:, 3), xbp_Wsq(:, 3)];

% Plot surface
surf(ax, G_surf, Q_surf, I_surf, ...
     EdgeColor='none', FaceColor=special_colour, FaceAlpha=0.25);

% Plot specific lines for theta_perturb = 0, 0.125, 0.25
% for idx = 1 : length(theta_old_SP)
for idx = 1 : 2
  % plotting vector
  % xplot = [xbp_gamma(idx, :); xbp_Wsq(idx, :)];
  xplot = [xbp_gamma_SP(idx, :); xbp_Wsq_SP(idx, :)];

  plot3(ax, xplot(:, 1), xplot(:, 2), xplot(:, 3), ...
        Color=special_colour, LineWidth=lw, LineStyle='-');
end

%--------------------%
%     Plot: DTCs     %
%--------------------%
% for idx = 1 : length(DTCs)
%   DTC_plot = DTCs{idx};
%   plot3(ax, DTC_plot(:, 1), DTC_plot(:, 2), DTC_plot(:, 3), ...
%         LineStyle='-', Color='k', LineWidth=lw);
% end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
% ax.XAxis.Limits = [0.0, 7.0];
% ax.YAxis.Limits = [0.0, 4.0];
% ax.ZAxis.Limits = [0.0, 20];
ax.XAxis.Limits = [0.8, 5.0];
ax.YAxis.Limits = [-1, 4];
ax.ZAxis.Limits = [0.0, 16.1];

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

%-----------------------%
%     Plot: Compass     %
%-----------------------%
% Add inset axis for orientation compax
ax_inset = axes(Parent=fig, Position=[0.6, 0.6, 0.2, 0.2]);

% Set data aspect ratio
daspect(ax_inset, [1 1 1])

% Coordinates for arrows
x_arrow = [[0, 1]; [0, 0]; [0, 0]];
y_arrow = [[0, 0]; [0, 1]; [0, 0]];
z_arrow = [[0, 0]; [0, 0]; [0, 1]];

% Hold axis
hold(ax_inset, 'on');

% Draw arrows
plot3(ax_inset, x_arrow(1, :), x_arrow(2, :), x_arrow(3, :), ...
      Color='k', LineWidth=ax.LineWidth);
plot3(ax_inset, y_arrow(1, :), y_arrow(2, :), y_arrow(3, :), ...
      Color='k', LineWidth=ax.LineWidth);
plot3(ax_inset, z_arrow(1, :), z_arrow(2, :), z_arrow(3, :), ...
      Color='k', LineWidth=ax.LineWidth);

% Turn off hold
hold(ax_inset, 'off');

% Turn off axis
box(ax_inset, 'off');
grid(ax_inset, 'off');
axis(ax_inset, 'off');

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'off');
grid(ax, 'off');
axis(ax, 'off');

% 3D plot view
% view_angle = [340, 8.5];
view_angle = [330, 4.5];
view(ax, view_angle);
view(ax_inset, view_angle);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig11a.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
