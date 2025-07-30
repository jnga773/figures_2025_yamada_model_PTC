% clear all; close all; clc;

%-------------------------------------------------------------------------%
%                         Read Periodic Orbit Data                        %
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig12_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Default colour order
colours = colororder();

%-------------------------------------------------------------------------%
%%                               Plot Data                               %%
%-------------------------------------------------------------------------%
% Default colour order (matplotlib)
colours = colororder();

% Setup figure
fig = figure(5); clf;
fig.Name = 'PTCs';
ax = gca();

% Axis dimensions
width = 3.0;
height = 5.0;

% Set figure size
set_figure_dimensions(width, height, scale=2);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%-----------------------------%
%     Plot: Patch and DTC     %
%-----------------------------%
% Fundamental domain
patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
    FaceAlpha=0.2, EdgeColor='none', ...
    HandleVisibility='off');

for idx = 1 : length(A_perturb)
  theta_perturb_plot = theta_perturb{idx};
  theta_new_plot     = theta_new{idx};

  label_str = sprintf('$A = %.4f$', A_perturb(idx));

  plot(ax, theta_perturb_plot, theta_new_plot, ...
       LineStyle='-', Color=colours(idx, :), ...
       DisplayName=label_str);
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = 0.0 : 0.25 : 0.5;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.125 : 0.5;

ax.YAxis.TickValues = -1.5 : 0.5 : 2.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -1.5 : 0.25 : 2.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\varphi_{\mathrm{d}}$');
% ylabel(ax, '$\vartheta_{\mathrm{n}}$');

% Turn off all tick labels
% ax.XAxis.TickLabels = {};
% ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [-1, 1];
ax.YAxis.Limits = [-1.5, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../pdf/fig12_DTCs.pdf';
% exportgraphics(fig, filename_out, ContentType='vector');
