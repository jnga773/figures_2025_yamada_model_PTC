%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig11_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Matplotlib colours
colours = {'#1f77b4';  % blue
           '#ff7f0e';  % orange
           '#2ca02c';  % green
           '#d62728';  % red
           '#9467bd';  % purple
           '#8c564b';  % brown
           '#e377c2';  % pink
           '#7f7f7f';  % gray
           '#bcbd22';  % yellow-green
           '#17becf'   % cyan
           };

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
width = 7.8;
height = 8.0;

% Set figure size
set_figure_dimensions(width, height, scale=1);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%-----------------------------%
%     Plot: Patch Below 0     %
%-----------------------------%
patch(ax, [-20, 20, 20, -20], [-20, -20, 0, 0], 'k', ...
      FaceAlpha=0.2, EdgeColor='none');

%-------------------------------------------%
%     Plot: Periodic Orbit and Manifold     %
%-------------------------------------------%
% Data points
gamma_plot = xbp_gamma_SP(3, :);
Wsq_plot   = xbp_Wsq_SP(3, :);

% Plot gamma_theta_old point
plot(ax, gamma_plot(1), gamma_plot(3), Marker='o', MarkerSize=5, ...
     MarkerEdgeColor='k', MarkerFaceColor=colours{3}, LineWidth=1.0);
% Plot Wsq point
plot(ax, Wsq_plot(1), Wsq_plot(3), Marker='o', MarkerSize=5, ...
     MarkerEdgeColor='k', MarkerFaceColor=colours{1}, LineWidth=1.0);

%--------------------%
%     Plot: DTCs     %
%--------------------%
% DTC colours
DTC_colours = {colours{9}, colours{4}, colours{5}};
% A_perturb values
DTC_A_perturb = [0.1, 0.724237, 10.0];
% plot theta
theta = 0 : 0.01 : 2 * pi;
circle = [cos(theta); sin(theta)];

% Plot
for idx = 2 : length(DTC_A_perturb)
  % Plot DTC
  DTC_plot = [gamma_plot(1); gamma_plot(3)] + (DTC_A_perturb(idx) * circle);
  plot(ax, DTC_plot(1, :), DTC_plot(2, :), LineWidth=2.0, LineStyle='-', ...
       Color=DTC_colours{idx});
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
daspect([1, 1, 1]);

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-10.0, 12.0];
ax.YAxis.Limits = [-2.0, 12.0];

%------------------------------%
%     Axis Ticks: Settings     %
%------------------------------%
ax.XAxis.TickDirection = 'in';
ax.XAxis.MinorTick = 'on';
ax.YAxis.TickDirection = 'in';
ax.YAxis.MinorTick = 'on';

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = -12 : 3 : 12;
ax.XAxis.MinorTickValues = -12 : 1 : 12;

% Y-Axis
ax.YAxis.TickValues = -12 : 3 : 12;
ax.YAxis.MinorTickValues = -12 : 1 : 12;

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
filename_out = './fig12a.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
