% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig11_data.mat', 'I_theta_SP');
load('../data_files/fig12_data.mat');

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

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Setup figure
fig = figure(4); clf;
fig.Name = 'DTCs';
ax = gca();

% Axis dimensions
width = 2.5;
height = 6;

% Set figure size
set_figure_dimensions(width, height, scale=1);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%---------------------%
%     Plot: Patch     %
%---------------------%
% % Fundamental domain
% patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
%     FaceAlpha=0.2, EdgeColor='none', ...
%     HandleVisibility='off');

%--------------------%
%     Plot: DTCs     %
%--------------------%
% Plot linewidth
lw = 1.5;
% Plot colours
DTC_colours = {colours{9}, colours{4}, colours{5}};

% MX line options
MX_colour = [0.0, 0.0, 0.0, 0.5];
MX_lw     = 1.0;
MX_ls     ='-';

for idx_A = 3
  % Read data
  x_plot = theta_perturb{idx_A};
  y_plot = theta_new{idx_A};

  % Plot highlight for max angle
  if idx_A == 2
    plot(ax, [0.125, 0.125], [-20, 20], ...
         Color=MX_colour, LineStyle=MX_ls, LineWidth=MX_lw);
  elseif idx_A == 3
    % Calculate angle which hits {I = 0}
    theta_die = acos(I_theta_SP(3) / A_perturb(idx_A));
    MX1 = ((1.5 * pi) + theta_die) / (2 * pi);
    MX2 = ((1.5 * pi) - theta_die) / (2 * pi);
    % % Plot lines
    % plot(ax, [MX1, MX1], [-20, 20], ...
    %      Color=MX_colour, LineStyle=MX_ls, LineWidth=MX_lw);
    % plot(ax, [MX2, MX2], [-20, 20], ...
    %      Color=MX_colour, LineStyle=MX_ls, LineWidth=MX_lw);
    % Plot patch
    patch(ax, [MX1, MX2, MX2, MX1], [-5, -5, 5, 5], ...
          'k', FaceAlpha=0.2, EdgeColor='none');
  end

  % Plot DTC
  for offset = -3 : 3
    plot(ax, x_plot-1, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
    plot(ax, x_plot, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
    plot(ax, x_plot+1, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
  end

end
%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = [0.0, 0.5, 1.0];
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = [0.25, 0.75];

ax.YAxis.TickValues = [-0.5, 0.0, 0.5, 1.0, 1.5];
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = [-0.25, 0.25, 0.75, 1.25];

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\varphi_{\mathrm{d}}$');
% ylabel(ax, '$\vartheta_{\mathrm{n}}$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1];
ax.YAxis.Limits = [-0.25, 1.25];

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
daspect(ax, [1, 1, 1]);

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = sprintf('./fig12_DTCs_%d.pdf', idx_A);
exportgraphics(fig, filename_out, ContentType='vector');
