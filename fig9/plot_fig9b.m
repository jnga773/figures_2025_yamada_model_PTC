% clear all; close all; clc;

%-------------------------------------------------------------------------%
%%                               Read Data                               %%
%-------------------------------------------------------------------------%
% Load data
load('../data_files/fig9_data.mat');

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
% plot_colours = {'#92b700';    % Green-Yellow
%                 '#e6b400';    % Yellow
%                 '#eb5e00';    % Orange
%                 '#d62728';    % Red
%                 '#e377c2';    % Pink
%                 '#bf42f5';    % Purple
%                 '#1f9ece'};   % Cyan
plot_colours = {'#92b700';    % Green-Yellow
                '#e6b400';    % Yellow
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5'};   % Purple

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
height = 7.0;

% Set figure size
set_figure_dimensions(width, height);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%----------------------------%
%     Plot: Other Things     %
%----------------------------%
% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), ...
      FaceAlpha=0.2, EdgeColor='none');

% Plot diagonal lines
plot(ax, [0, 1], [0, 1], LineStyle='-', LineWidth=1.5, ...
     Color=colours(3, :));

% Grey lines at theta_old = 0 and 0.3
xline(ax, 0.5148, LineStyle='-', LineWidth=1, ...
      Color=[0, 0, 0, 0.5]);

%--------------------%
%     Plot: PTCs     %
%--------------------%
% Linewidth
lw = 1.5;

% Plot all PTCs
for idx = 1 : length(A_perturb)
  fprintf('A_p = %.3f\n', A_perturb(idx));

  if idx == 2
    temp = -1;
  else
    temp = 0;
  end

  % Plot
  plot(ax, theta_old{idx}, theta_new{idx}+temp, ...
       LineWidth=lw, LineStyle='-', ...
       Color=plot_colours{idx});
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------------%
%     Data Aspect Ratio     %
%---------------------------%
ax.DataAspectRatio = [1, 1, 1];

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = -0.5 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = -0.5 : 0.25 : 1.0;

ax.YAxis.TickValues = -0.5 : 0.5 : 2.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.5 : 0.25 : 2.5;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% Axis labels
% xlabel(ax, '$\theta_{\mathrm{o}}$');
% ylabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all tick labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1.0];
ax.YAxis.Limits = [-0.25, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../pdf/fig9b_I_PTCs.pdf';
exportgraphics(fig, filename_out, ContentType='vector');
