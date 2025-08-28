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
fig = figure(3); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
width = 3.5;
height = 2.0;

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
plot(ax, xbp_PO1(:, 1), xbp_PO1(:, 2), ...
     Color=colour_PO_transparent, LineWidth=lw);
plot(ax, xbp_PO2(:, 1), xbp_PO2(:, 2), ...
     Color=colour_PO_transparent, LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot(ax, xbp_gamma(:, 1), xbp_gamma(:, 2), ...
     Color=colour_PO, LineWidth=lw, LineStyle='-');

%-------------------------------%
%     Plot: Stable Manifold     %
%-------------------------------%
% Plot original stable manifold
plot(ax, Wqs1(:, 1), Wqs1(:, 2), ...
     Color=colour_Wsq_transparent, LineWidth=lw);
plot(ax, Wqs2(:, 1), Wqs2(:, 2), ...
     Color=colour_Wsq_transparent, LineWidth=lw);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 2), ...
     Color=colour_Wsq, LineWidth=lw, LineStyle='-');

%---------------------------------%
%     Plot: Equilibrium Point     %
%---------------------------------%
% Plot equilibrium point
plot(ax, xpos(1), xpos(2), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Plot specific lines for theta_perturb = 0, 0.125, 0.25
for idx = [1, 2]
  % plotting vector
  xplot = [xbp_gamma_SP(idx, :); xbp_Wsq_SP(idx, :)];

  plot(ax, xplot(:, 1), xplot(:, 2), ...
       Color=colour_special, LineWidth=lw, LineStyle='-');
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = [];

% Z-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = [];

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.75, 3.0];
ax.YAxis.Limits = [0.0, 2.0];

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$G$');
% ylabel(ax, '$I$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'off');

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = './fig11a_inset2.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
