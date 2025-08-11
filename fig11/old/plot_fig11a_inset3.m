%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s', 'xpos', 'tbp_PO', 'T_PO');
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
% Special colour
special_colour = colours(7, :);

% Plot original periodic orbit
plot(ax, xbp_PO1(:, 1), xbp_PO1(:, 3), ...
     Color=colours(3, :), ...
     LineWidth=lw);
plot(ax, xbp_PO2(:, 1), xbp_PO2(:, 3), ...
     Color=colours(3, :), ...
     LineWidth=lw);

% Plot special section
plot(ax, xbp_gamma(:, 1), xbp_gamma(:, 3), ...
     Color=special_colour, LineWidth=lw, LineStyle='-');

%-------------------------------%
%     Plot: Stable Manifold     %
%-------------------------------%
% Plot original stable manifold
plot(ax, Wqs1(:, 1), Wqs1(:, 3), ...
     Color=colours(1, :), ...
     LineWidth=lw);
plot(ax, Wqs2(:, 1), Wqs2(:, 3), ...
     Color=colours(1, :), ...
     LineWidth=lw);

% Plot special section
plot(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 3), ...
     Color=special_colour, LineWidth=lw, LineStyle='-');

%---------------------------------%
%     Plot: Equilibrium Point     %
%---------------------------------%
% Plot equilibrium point
plot(ax, xpos(1), xpos(3), ...
     Marker='o', MarkerSize=4, ...
     MarkerFaceColor='r', MarkerEdgecolor='k', LineWidth=0.25);

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Plot specific lines for theta_perturb = 0, 0.125, 0.25
for idx = [1, 2]
  % plotting vector
  xplot = [xbp_gamma_SP(idx, :); xbp_Wsq_SP(idx, :)];

  plot(ax, xplot(:, 1), xplot(:, 3), ...
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
ax.YAxis.Limits = [0.0, 5.0];

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
%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig11a_inset3.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
