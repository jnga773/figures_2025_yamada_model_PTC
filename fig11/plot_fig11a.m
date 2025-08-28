%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig11_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO = '#2ca02c';
% Stable manifold colour
colour_Wsq = '#1f77b4';

% Transparent versions
colour_PO_transparent  = [hex2rgb(colour_PO), 0.2];
colour_Wsq_transparent = [hex2rgb(colour_Wsq), 0.2];

% Highlight colour
colour_special = '#e377c2';

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

% Plot original periodic orbit
plot3(ax, xbp_PO1(:, 1), xbp_PO1(:, 2), xbp_PO1(:, 3), ...
      Color=colour_PO_transparent, LineWidth=lw);
plot3(ax, xbp_PO2(:, 1), xbp_PO2(:, 2), xbp_PO2(:, 3), ...
      Color=colour_PO_transparent, LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_gamma(:, 1), xbp_gamma(:, 2), xbp_gamma(:, 3), ...
      Color=colour_PO, LineWidth=lw, LineStyle='-');

%-------------------------------%
%     Plot: Stable Manifold     %
%-------------------------------%
% Plot original stable manifold
plot3(ax, Wqs1(:, 1), Wqs1(:, 2), Wqs1(:, 3), ...
      Color=colour_Wsq_transparent, LineWidth=lw);
plot3(ax, Wqs2(:, 1), Wqs2(:, 2), Wqs2(:, 3), ...
      Color=colour_Wsq_transparent, LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 2), xbp_Wsq(:, 3), ...
      Color=colour_Wsq, LineWidth=lw, LineStyle='-');

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
     EdgeColor='none', FaceColor=colour_special, FaceAlpha=0.25);

% Plot specific lines for theta_perturb = 0, 0.125, 0.25
% for idx = 1 : length(theta_old_SP)
for idx = 1 : 2
  % plotting vector
  % xplot = [xbp_gamma(idx, :); xbp_Wsq(idx, :)];
  xplot = [xbp_gamma_SP(idx, :); xbp_Wsq_SP(idx, :)];

  plot3(ax, xplot(:, 1), xplot(:, 2), xplot(:, 3), ...
        Color=colour_special, LineWidth=lw, LineStyle='-');
end

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
view_angle = [330, 4.5];
view(ax, view_angle);
view(ax_inset, ax.View);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = './fig11a.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
