close all; clear; clc;

%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s', 'xpos');
load('../data_files/fig11_data.mat', 'theta_old', 'A_perturb', 'theta_perturb');

%----------------------%
%     Plot Colours     %
%----------------------%
% Default colour order
colours = colororder();

%%
%-------------------------------------------------------------------------%
%                         Plot Data: 2D Comparison                        %
%-------------------------------------------------------------------------%
colours = colororder();

% Setup figure
fig = figure(2); clf;
fig.Name = 'Intersection with Ws(q) and Perturbed Gamma';
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

%--------------------------------------------%
%     Plot: Highlight theta_perturb Area     %
%--------------------------------------------%
% Highlight area over 0.0 <= theta_perturb <= 0.5 pi
patch(ax(1), [0, 0, 0.25, 0.25], [0, 1, 1, 0], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

%------------------------------------------%
%     Plot: theta_perturb vs theta_old     %
%------------------------------------------%
% Left axis
yyaxis(ax, 'left');

% Plot: theta_old vs theta_perturb
plot(ax(1), theta_perturb, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
plot(ax(1), theta_perturb, theta_old-1, Color=colours(1, :), LineStyle='-', LineWidth=2.5);
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

%----------------------------%
%     Axis Ticks: Axis 1     %
%----------------------------%
% X-Axis: theta_perturb
ax.XAxis.MinorTick = 'on';
ax.XAxis.TickValues = -0.25 : 0.25 : 1.5;
ax.XAxis.MinorTickValues = -0.25 : 0.125 : 1.5;

% Y-Axis 1: theta_old
ax.YAxis(1).MinorTick = 'on';
ax.YAxis(1).TickValues = 0.0 : 0.25 : 1.0;
ax.YAxis(1).MinorTickValues = 0.0 : 0.125 : 1.0;

% Y-Axis 2: A_perturb
ax.YAxis(2).MinorTick = 'on';
ax.YAxis(2).TickValues = 0.0 : 3.0 : 15.0;
ax.YAxis(2).MinorTickValues = 0.0 : 1.0 : 15.0;

%---------------------%
%     Tick Labels     %
%---------------------%
ax.XAxis.TickLabels = {};
ax.YAxis(1).TickLabels = {};
ax.YAxis(2).TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
% X-Axis: theta_perturb
ax.XAxis.Limits = [-0.125, 1.0];

% Y-Axis 1: theta_old
ax.YAxis(1).Limits = [0.0, 1.0];

% Y-Axis 2: A_perturb
ax.YAxis(2).Limits = [0.0, 15.0];

%---------------------%
%     Axis Labels     %
%---------------------%
% % X-Axis: theta_perturb
% ax.XAxis.Label.String = '$\varphi_{\mathrm{p}}$';
% 
% % Y-Axis 1: theta_perturb
% ax.YAxis(1).Label.String = '$\vartheta_{\mathrm{o}}$';
% 
% % Y-Axis 2: A_perturb
% ax.YAxis(2).Label.String = '$A$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

% Save figure
filename_out = '../pdf/fig11b.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
