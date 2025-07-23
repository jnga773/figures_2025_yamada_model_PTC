%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
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
fig = figure(1); clf;
fig.Name = 'Periodic Orbit Phase Portrait (3D)';
ax = gca();

% Axis dimensions
width = 7.5;
height = 5.0;

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

% Plot highlighted sections along \Gamma and W^{s}(q)
plot3(ax, xbp_gamma(:, 1), xbp_gamma(:, 2), xbp_gamma(:, 3), ...
      Color=colour, LineWidth=lw, LineStyle='-');
plot3(ax, xbp_Wsq(:, 1), xbp_Wsq(:, 2), xbp_Wsq(:, 3), ...
      Color=colour, LineWidth=lw, LineStyle='-');

% Create surface plots
G_surf = [xbp_gamma(:, 1), xbp_Wsq(:, 1)];
Q_surf = [xbp_gamma(:, 2), xbp_Wsq(:, 2)];
I_surf = [xbp_gamma(:, 3), xbp_Wsq(:, 3)];

% Plot surface
surf(ax, G_surf, Q_surf, I_surf, ...
     EdgeColor='none', FaceColor=colour, FaceAlpha=0.5);

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
% ax.XAxis.Limits = [0.0, 7.0];
% ax.YAxis.Limits = [0.0, 4.0];
% ax.ZAxis.Limits = [0.0, 20];
ax.XAxis.Limits = [0.85, 5.0];
ax.YAxis.Limits = [0.2, 4.0];
ax.ZAxis.Limits = [0.0, 16.1];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickDirection = 'in';
ax.XAxis.TickValues = 0.0 : 2.0 : 8.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 1.0 : 7.0;

% Y-Axis
ax.YAxis.TickDirection = 'in';
ax.YAxis.TickValues = 0.0 : 2.0 : 4.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 1.0 : 4.0;

% Z-Axis
ax.ZAxis.TickDirection = 'in';
ax.ZAxis.TickValues = 0.0 : 4.0 : 20.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 0.0 : 2 : 20.0;

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
grid(ax, 'on');

% Grid lines
ax.GridLineWidth = 0.5;
ax.GridColor = 'black';
ax.GridAlpha = 0.25;

% 3D plot view
% view(45, 6.0);
% view(-45, 6);
view(-3.5, 4);

%---------------------%
%     Save Figure     %
%---------------------%
filename_out = '../pdf/fig11a.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
