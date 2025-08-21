clear all; close all; clc;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     % . 
%-------------------%
% Load big PTC scan data
load('../data_files/fig7_data.mat');

% Shift data_hole_lt1 start data point
data_hole_lt1 = fix_gap(data_before_hole, data_hole_lt1);

%--------------------------------%
%     Coordinates for 'Hole'     %
%--------------------------------%
% G-direction
intersection.theta_old = 0.326074664455562;
intersection.A_perturb = 0.543525459141593;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
plot_colours = {'#92b700';    % Green-Yellow
                '#e6b400';    % Yellow
                '#eb5e00';    % Orange
                '#d62728';    % Red
                '#e377c2';    % Pink
                '#bf42f5';    % Purple
                '#1f9ece'};   % Cyan

%-----------------------------------%
%     Sort Out Single Plot Data     %
%-----------------------------------%
A_perturb_specific = [0.05, 0.1, 0.15, 0.5435, 1.0, 1.5, 2.0];

% Find specific curves
[theta_old_plot, A_perturb_plot, theta_new_plot] = ...
    find_PTC_curves('../data_files/fig7_data.mat', A_perturb_specific);

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(1); clf;
fig.Name = 'PTC Scans: Gain';
ax  = axes;

% Axis dimensions
width = 7.8;
height = 6.0;

% Set figure size
set_figure_dimensions(width, height, scale=3);

% Set axis linewidth
ax.LineWidth = 0.8;

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'on');

%-------------------------------------------------%
%     Plot: Stable Manifold Intersection Pole     %
%-------------------------------------------------%
plot3(ax, [intersection.theta_old, intersection.theta_old], ...
     [intersection.A_perturb, intersection.A_perturb], ...
     [-5, 5], ...
     Color='k', LineWidth=2.5, LineStyle='-');

%--------------------------%
%     Surface Settings     %
%--------------------------%
% Set colour map
cmap = colormap(scale_colour_map(2.0));

% Shading of surface
shading(ax, 'interp');

% % Lighting
% light(ax, Style='infinite', Position=[intersection.theta_old, intersection.A_perturb, 20]);

% Facealpha
% facealpha = 0.75;
facealpha = 1.0;

%---------------------------%
%     Plot: Before Hole     %
%---------------------------%
% Surface: Before hole
[X, Y, Z] = pad_data(data_before_hole);
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z-1, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

%--------------------------%
%     Plot: After Hole     %
%--------------------------%
% Surface: After hole
[X, Y, Z] = pad_data(data_after_hole);
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z+1, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

%------------------------%
%     Plot: Hole LHS     %
%------------------------%
% Surface: Hole (theta_old > 1)
[X, Y, Z] = pad_data(data_hole_gt1, 'gt1');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

% Plot rim
[theta_old_rim, A_perturb_rim, theta_new_rim] = find_rim_data(X, Y, Z, 'gt1');
plot3(ax, theta_old_rim, A_perturb_rim, theta_new_rim, ...
      LineStyle='-', Color='k', LineWidth=2.0);

%------------------------%
%     Plot: Hole RHS     %
%------------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = pad_data(data_hole_lt1, 'lt1');
surf(ax, X, Y, Z, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');
surf(ax, X, Y, Z+1.0, EdgeColor='none', FaceColor='interp', MeshStyle='row', ...
     FaceAlpha=facealpha, FaceLighting='flat');

% Plot rim
[theta_old_rim, A_perturb_rim, theta_new_rim] = find_rim_data(X, Y, Z, 'lt1');
plot3(ax, theta_old_rim, A_perturb_rim, theta_new_rim, ...
      LineStyle='-', Color='k', LineWidth=2.0);
plot3(ax, theta_old_rim, A_perturb_rim, theta_new_rim+1, ...
      LineStyle='-', Color='k', LineWidth=2.0);

%--------------------------%
%     Plot: PTC Curves     %
%--------------------------%
% Linewidth
lw = 2.0;
% 
% Plot all PTCs
for idx= 1 : length(A_perturb_specific)
  % Plot specific PTCs over the lower surface
  plot3(ax, theta_old_plot{1, idx}, A_perturb_plot{1, idx}, theta_new_plot{1, idx}, ...
        LineWidth=lw, LineStyle='-', Color=plot_colours{idx});
  % Plot other side of hole
  if idx == 4 || idx == 5
    plot3(ax, theta_old_plot{2, idx}, A_perturb_plot{2, idx}, theta_new_plot{2, idx}, ...
          LineWidth=lw, LineStyle='-', Color=plot_colours{idx});
  end

  % Plot specific PTCs over the upper surface
  if idx <= 3
    plot3(ax, theta_old_plot{1, idx}, A_perturb_plot{1, idx}, theta_new_plot{1, idx}-1, ...
          LineWidth=lw, LineStyle='-', Color=plot_colours{idx});
  else
    plot3(ax, theta_old_plot{1, idx}, A_perturb_plot{1, idx}, theta_new_plot{1, idx}+1, ...
          LineWidth=lw, LineStyle='-', Color=plot_colours{idx});
  end
end

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 2.0];
ax.ZAxis.Limits = [0.0, 3.25];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 1.0 : 2.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 0.5 : 2.0;

% Z-Axis
ax.ZAxis.TickValues = 0.0 : 0.5 : 3.0;
ax.ZAxis.MinorTick = 'on';
ax.ZAxis.MinorTickValues = 0.0 : 0.25 : 3.0;

%------------------------------%
%     Axis and Tick Labels     %
%------------------------------%
% % Axis labels
% xlabel(ax, '$\theta_{\mathrm{o}}$');
% ylabel(ax, '$A_{\mathrm{p}}$');
% zlabel(ax, '$\theta_{\mathrm{n}}$');

% Turn off all axis labels
ax.XAxis.TickLabels = {};
ax.YAxis.TickLabels = {};
ax.ZAxis.TickLabels = {};

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
grid(ax, 'on');

% axis(ax, 'off');

%---------------------%
%     Save Figure     %
%---------------------%
view(315, 15);

filename_out = './fig7_G_PTC_full.png';
% exportgraphics(fig, filename_out, ContentType='image', Resolution=1000);

%-------------------------------------------------------------------------%
%%                               FUNCTION                                %%
%-------------------------------------------------------------------------%
function colourmap_out = scale_colour_map(scale_factor)
  % colourmap_out = scale_colour_map(scale_factor)
  % 
  % Rescale the colour map.

  n_colours_in = 2048;
  n_colours_out = 1024;

  % Get colour map
  colour_map = parula(n_colours_in);

  % Create a linear ramp the size of the colormap we actually want
  t = linspace(0,1,n_colours_out);
  % Apply whatever transform you like to the ramp
  t2 = t .^ scale_factor;

  % Use that to scale the big linear colormap into the small stretched one.
  colourmap_out = colour_map(1+floor((n_colours_in-1)*t2'),:);
end
