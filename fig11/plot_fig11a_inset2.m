%%
%-------------------------------------------------------------------------%
%                              Calculate Data                             %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s', 'xpos');
load('../data_files/fig11_data.mat', 'xbp_gamma', 'xbp_Wsq');

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
fig = figure(3); clf;
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
lw = 2.5;
colour = colours(7, :);

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_gamma(:, 1), xbp_gamma(:, 2), xbp_gamma(:, 3), ...
      Color=colour, LineWidth=lw, LineStyle='-');
plot3(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 2), xbp_Wsq(:, 3), ...
      Color=colour, LineWidth=lw, LineStyle='-');

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
ax.YAxis.Limits = [0.0, 2.0];
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

% 3D plot view
% view(45, 6.0);
% view(-41, -4);
view(-70, 6);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig11a_inset2.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
