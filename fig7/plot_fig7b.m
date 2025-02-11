% Clear stuff
clear all; close all;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Load big PTC scan data
load('../data_files/fig7_data.mat');

%--------------------------------%
%     Coordinates for 'Hole'     %
%--------------------------------%
% G-direction
intersection.theta_old = 0.3176;
intersection.A_perturb = 0.5576;

%----------------------%
%     Plot Colours     %
%----------------------%
% Plot colours
% Green     (#2ca02c) = [ 44, 160,  44] ./ 255
% Chartreus (#bcbd22) = [188, 189,  34] ./ 255
% Yellow    (#fafa2a) = [250, 250,  42] ./ 255
% Orange    (#ff7f0e) = [255, 126,  14] ./ 255
% Red       (#d62728) = [214,  39,  40] ./ 255
% Pink      (#e38ab7) = [227, 138, 183] ./ 255
% Purple    (#9467bd) = [148, 103, 189] ./ 255
% Cyan      (#1bc3cc) = [ 27, 195, 204 ./ 255

% Plot colours
plot_colours = {[188, 189,  34] ./ 255;
                [250, 250,  42] ./ 255;
                [255, 126,  14] ./ 255;
                [214,  39,  40] ./ 255;
                [227, 138, 183] ./ 255;
                [148, 103, 189] ./ 255;
                [ 27, 195, 204] ./ 255};

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Colourmap     %
%-------------------%
% Plotting colours
colours = colororder();

%-------------------------%
%     Figure Settings     %
%-------------------------%
fig = figure(7); clf;
fig.Name = 'PTC Scans: Gain';
ax = gca();

% Axis dimensions
width = 8.0;
height = 6.4;

% Add set_figure_dimensions() function to path
% addpath('../');

% Set figure size
set_figure_dimensions(width, height);

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
%     Plot: PTC Curves     %
%--------------------------%
% Linewidth
lw = 1.5;

plot_A_perturb = [0.05, 0.1, 0.15, 0.55, 1.0, 1.5, 2.0];

% Find plotting indices
plot_idx = zeros(length(plot_A_perturb), 1);
for i = 1 : length(plot_A_perturb)
  plot_idx(i) = find(round(A_perturb, 3) == plot_A_perturb(i));
end

% Plot all PTCs
for i = 1 : length(plot_idx)
  idx = plot_idx(i);

  fprintf('A_p = %.3f\n', A_perturb(idx));

  % A_perturb plot data
  A_plot_lt1 = A_perturb(idx) * ones(length(theta_old_lt1{idx}));
  A_plot_gt1 = A_perturb(idx) * ones(length(theta_old_gt1{idx}));

  % Plot
  plot3(ax, theta_old_lt1{idx}, A_plot_lt1, theta_new_lt1{idx}, Color=plot_colours{i}, LineStyle='-');
  plot3(ax, theta_old_gt1{idx}-1, A_plot_gt1, theta_new_gt1{idx}, Color=plot_colours{i}, LineStyle='-');
end

%--------------------------%
%     Surface Settings     %
%--------------------------%
% Set colour map
colormap(scale_colour_map(0.9));

% Shading of surface
shading(ax, 'interp');

%-----------------------%
%     Plot: Surface     %
%-----------------------%
% Surface: Hole (theta_old < 1)
[X, Y, Z] = pad_data(data_hole_lt1, 0, 'lt1');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: Hole (theta_old > 1)
[X, Y, Z] = pad_data(data_hole_gt1, 0, 'gt1');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: Before hole
[X, Y, Z] = pad_data(data_before_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

% Surface: After hole
[X, Y, Z] = pad_data(data_after_hole, 0, 'none');
surf(ax, X, Y, Z, EdgeColor='interp', FaceColor='interp', MeshStyle='row');

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0.0, 1.0];
ax.YAxis.Limits = [0.0, 2.0];
ax.ZAxis.Limits = [0.0, 2.5];

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

%---------------------%
%     Save Figure     %
%---------------------%
% View point
view(135, 15);

% Filename
filename_out = '../pdf/fig7b_G_PTC_surface_2.png';
% exportgraphics(fig, filename_out, ContentType='imagwe', Resolution=750);
