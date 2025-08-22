%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig11_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO   = '#2ca02c';
% Stable manifold colour
colour_Wsq  = '#1f77b4';

% Plot colours
DTC_colours = {'#bcbd22';
               '#d62728';
               '#9467bd'};

%--------------------------%
%     Calculate Things     %
%--------------------------%
% Data points
gamma_plot = xbp_gamma_SP(3, :);
Wsq_plot   = xbp_Wsq_SP(3, :);

% A_perturb values
DTC_A_perturb = [0.1, 0.724237, 10.0];

% Create circle data to plot DTC range
theta = 0.0 : 0.001 : 1.0;
circle = [cos(theta * (2 * pi)); sin(theta * (2 * pi))];

% Create DTC circles
DTC_circles = cell(1, length(DTC_A_perturb));
for idx = 1 : length(DTC_A_perturb)
  DTC_circles{idx} = [gamma_plot(1); gamma_plot(3)] + DTC_A_perturb(idx) * circle;
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
width = 3.0;
height = 3.0;

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

%--------------------%
%     Plot: DTCs     %
%--------------------%
for idx = 1 : 2
  DTC_plot = DTC_circles{idx};

  % Plot DTC
  plot(ax, DTC_plot(1, :), DTC_plot(2, :), LineWidth=2.0, LineStyle='-', ...
       Color=DTC_colours{idx});
end

%-------------------------------------------%
%     Plot: Periodic Orbit and Manifold     %
%-------------------------------------------%
% Data points
gamma_plot = xbp_gamma_SP(3, :);
Wsq_plot   = xbp_Wsq_SP(3, :);

% Plot gamma_theta_old point
plot(ax, gamma_plot(1), gamma_plot(3), Marker='o', MarkerSize=5, ...
     MarkerEdgeColor='k', MarkerFaceColor=colours_PO, LineWidth=1.0);
% Plot Wsq point
plot(ax, Wsq_plot(1), Wsq_plot(3), Marker='o', MarkerSize=5, ...
     MarkerEdgeColor='k', MarkerFaceColor=colours_Wsq, LineWidth=1.0);

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
% ax.XAxis.Limits = [0.4, 2.0];
% ax.YAxis.Limits = [0.6, 2.4];

ax.XAxis.Limits = [gamma_plot(1) - DTC_A_perturb(2) - 0.05, ...
                   gamma_plot(1) + DTC_A_perturb(2) + 0.05];
ax.YAxis.Limits = [gamma_plot(3) - DTC_A_perturb(2) - 0.05, ...
                   gamma_plot(3) + DTC_A_perturb(2) + 0.05];

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
ax.XAxis.TickValues = [];
ax.XAxis.MinorTickValues = [];

% Y-Axis
ax.YAxis.TickValues = [];
ax.YAxis.MinorTickValues = [];

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
filename_out = './fig12a_inset.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
