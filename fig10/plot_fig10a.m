% Clear stuff
clear all; close all;

%-------------------------------------------------------------------------%
%%                             Read Data                                 %%
%-------------------------------------------------------------------------%
%-------------------%
%     Read Data     %
%-------------------%
% Load big PTC scan data
load('../data_files/fig10_data.mat');

%--------------------------------%
%     Coordinates for 'Hole'     %
%--------------------------------%
% G-direction
intersection.theta_old = 0.5148;
intersection.A_perturb = 4.0737;

%-------------------------------------------------------------------------%
%%                             Plot Data                                 %%
%-------------------------------------------------------------------------%
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
%     Surface Settings     %
%--------------------------%
% Set colour map
colormap(scale_colour_map(2.0));

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
for i = 1 : length(data_before_hole.theta_new)
  data_before_hole.theta_new{i} = data_before_hole.theta_new{i} - 1.0;
end
data_before_hole.theta_old
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
ax.YAxis.Limits = [0.0, 25.0];
ax.ZAxis.Limits = [0.25, 3.0];

%--------------------%
%     Axis Ticks     %
%--------------------%
% X-Axis
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 1.0;

% Y-Axis
ax.YAxis.TickValues = 0.0 : 5.0 : 25.0;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = 0.0 : 2.5 : 25.0;

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
view(315, 15);

filename_out = '../pdf/fig10a_I_PTC_surface_1.png';
% exportgraphics(fig, filename_out, ContentType='image', Resolution=1000);
