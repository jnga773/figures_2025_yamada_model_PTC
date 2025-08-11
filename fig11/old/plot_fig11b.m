%%
%-------------------------------------------------------------------------%
%                                Read Data                                %
%-------------------------------------------------------------------------%
load('../data_files/fig2_data.mat', 'xbp_PO', 'Wq_s', 'xpos');
load('../data_files/fig11_data.mat', 'theta_old', 'A_perturb', 'theta_perturb');

%-------------------%
%     Sort Data     %
%-------------------%
% Sort indices
[~, sort_idx] = sort(theta_old);
theta_old     = theta_old(sort_idx);
theta_perturb = theta_perturb(sort_idx);

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
fig = figure(4); clf;
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
patch(ax, [0, 0, 0.25, 0.25], [0, 1, 1, 0], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

%------------------------------------------%
%     Plot: theta_perturb vs theta_old     %
%------------------------------------------%
% Linewidth
lw = 1.0;

% Plot: theta_old vs theta_perturb
plot(ax, theta_perturb, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=lw);
plot(ax, theta_perturb, theta_old-1, Color=colours(1, :), LineStyle='-', LineWidth=lw);
plot(ax, theta_perturb+2, theta_old, Color=colours(1, :), LineStyle='-', LineWidth=lw);

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

% Y-Axis: theta_old
ax.YAxis.MinorTick = 'on';
ax.YAxis.TickValues = 0.0 : 0.25 : 1.0;
ax.YAxis.MinorTickValues = 0.0 : 0.125 : 1.0;

%---------------------%
%     Tick Labels     %
%---------------------%
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
% X-Axis: theta_perturb
ax.XAxis.Limits = [-0.125, 1.0];

% Y-Axis: theta_old
ax.YAxis.Limits = [0.0, 1.0];

%---------------------%
%     Axis Labels     %
%---------------------%
% % X-Axis: theta_perturb
% ax.XAxis.Label.String = '$\varphi_{\mathrm{p}}$';
% 
% % Y-Axis: theta_perturb
% ax.YAxis.Label.String = '$\vartheta_{\mathrm{o}}$';

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');

% Save figure
filename_out = '../pdf/fig11b.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
