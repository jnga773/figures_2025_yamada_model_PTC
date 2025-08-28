%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s');
load('../data_files/fig11_data.mat', 'xbp_gamma_SP', 'xbp_Wsq_SP');

%----------------------%
%     Plot Colours     %
%----------------------%
% Periodic orbit colour
colour_PO  = '#2ca02c';
% Stable manifold colour
colour_Wsq = '#1f77b4';

% Transparent versions
colour_PO_transparent  = [hex2rgb(colour_PO), 0.4];
colour_Wsq_transparent = [hex2rgb(colour_Wsq), 0.4];

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
% Plot gamma_theta_old point
plot(ax, gamma_plot(1), gamma_plot(3), Marker='o', MarkerSize=4, ...
     MarkerEdgeColor='k', MarkerFaceColor=colour_PO, LineWidth=0.8);
% Plot Wsq point
plot(ax, Wsq_plot(1), Wsq_plot(3), Marker='o', MarkerSize=4, ...
     MarkerEdgeColor='k', MarkerFaceColor=colour_Wsq, LineWidth=0.8);

% % Plot 2d periodic orbit data
% plot(ax, xbp_PO(:, 1), xbp_PO(:, 3), ...
%      Color=colour_PO_transparent, LineStyle='-');
% plot(ax, Wq_s(:, 1), Wq_s(:, 3), ...
%      Color=colour_Wsq_transparent, LineStyle='-');

%--------------------%
%     Plot: DTCs     %
%--------------------%
for idx = 1 : length(DTC_A_perturb)
  DTC_plot = DTC_circles{idx};

  if idx == 3
    % Grey out section for A_perturb = 10
    idx1 = DTC_plot(2, :) < 0.0;

    % Plot DTC
    plot(ax, DTC_plot(1, idx1), DTC_plot(2, idx1), LineWidth=2.0, LineStyle=':', ...
         Color=DTC_colours{idx});

    % Regular plot otherwise
    idx2 = DTC_plot(2, :) >= 0.0;
    % Find where there is a big difference
    DTC_plot = DTC_plot(:, idx2);
    [~, max_idx] = max(diff(DTC_plot(1, :)));

    % Plot DTC
    plot(ax, DTC_plot(1, 1:max_idx), DTC_plot(2, 1:max_idx), LineWidth=2.0, LineStyle='-', ...
         Color=DTC_colours{idx});
    plot(ax, DTC_plot(1, max_idx+1:end), DTC_plot(2, max_idx+1:end), LineWidth=2.0, LineStyle='-', ...
         Color=DTC_colours{idx});

  else
    % Plot DTC
    plot(ax, DTC_plot(1, :), DTC_plot(2, :), LineWidth=2.0, LineStyle='-', ...
         Color=DTC_colours{idx});
  end
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
ax.XAxis.MinorTickValues = -12 : 1.5 : 12;

% Y-Axis
ax.YAxis.TickValues = -12 : 3 : 12;
ax.YAxis.MinorTickValues = -12 : 1.5 : 12;

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
