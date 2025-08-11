% clear all; close all;

%% READ DATA
run_in = 'run08_PR_DTC_scan';

%------------------------------%
%     Read Sub-Directories     %
%------------------------------%
% Folder name
dirs_P = sprintf('./data/%s/', run_in);
% List all directories
dirs_P = dir(dirs_P);
% Remove ./ and ../
dirs_P = dirs_P(~ismember({dirs_P.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dirs_P = {dirs_P.name};

%------------------------------------%
%     Cycle Through and Read DTC     %
%------------------------------------%
% Empty cells for data
A_perturb_data     = zeros(length(dirs_P), 4);
theta_old_data     = zeros(length(dirs_P), 4);
theta_new_data     = cell(length(dirs_P), 4);
theta_perturb_data = cell(length(dirs_P), 4);

% Cycle through A_perturb sub-directories
for idx_P = 1 : length(dirs_P)
  % Run name
  P_sub_run = {run_in, dirs_P{idx_P}};

  % Read folders inside dirs_P
  dirs_DTC = dir(sprintf('./data/%s/%s/', run_in, dirs_P{idx_P}));
  dirs_DTC = dirs_DTC(~ismember({dirs_DTC.name}, {'.', '..', '.DS_Store'}));
  dirs_DTC = {dirs_DTC.name};

  % Cycle through DTC sub-sub directories and read data
  for idx_DTC = 1 : length(dirs_DTC)
    DTC_sub_run = {run_in, dirs_P{idx_P}, dirs_DTC{idx_DTC}};

    % Read bifurcation data
    bd_read = coco_bd_read(DTC_sub_run);

    % Get values
    A_perturb_read     = coco_bd_val(bd_read, 1, 'A_perturb');
    theta_old_read     = coco_bd_val(bd_read, 1, 'theta_old');
    theta_new_read     = coco_bd_col(bd_read, 'theta_new');
    theta_perturb_read = coco_bd_col(bd_read, 'theta_perturb');
    
    % Update arrays
    A_perturb_data(idx_P, idx_DTC)     = A_perturb_read;
    theta_old_data(idx_P, idx_DTC)     = theta_old_read;
    theta_new_data{idx_P, idx_DTC}     = theta_new_read;
    theta_perturb_data{idx_P, idx_DTC} = theta_perturb_read;
  end

end

%% COMBINE DATA
idx_A = 1;

[x_plot, y_plot] = sort_mod_data(theta_perturb_data, theta_new_data, idx_A);

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
width = 6;
height = 8.5;

% Set figure size
set_figure_dimensions(width, height, scale=1);

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

lw = 1.5;

DTC_colours = {colours(9, :); colours(4, :); colours(5, :)};
for idx_A = 1 : 3
  % Read data
  [x_plot, y_plot] = sort_mod_data(theta_perturb_data, theta_new_data, idx_A);

  % Plot DTC
  for offset = -1 : 1
    plot(ax, x_plot, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
    plot(ax, x_plot+1, y_plot+offset, LineStyle='-', Color=DTC_colours{idx_A}, LineWidth=lw);
  end
end

hold(ax, 'off');

daspect(ax, [1, 1, 1]);

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%--------------------%
%     Axis Ticks     %
%--------------------%
ax.XAxis.TickValues = 0.0 : 0.5 : 1.0;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 0.0 : 0.25 : 0.5;

ax.YAxis.TickValues = -0.5 : 0.5 : 1.5;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = -0.25 : 0.25 : 1.25;

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

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');

%----------------------%
%      Save Figure     %
%----------------------%
% Filename
filename_out = '../pdf/fig12_DTCs.pdf';
exportgraphics(fig, filename_out, ContentType='vector');

%%
function [x_out, y_out] = sort_mod_data(theta_perturb_data_in, theta_new_data_in, idx_A_in)
  % Cycle through theta_perturb directories and append data
  theta_perturb_plot = [];
  theta_new_plot     = [];

  % Get dimensions of data array
  dir_P_size = size(theta_perturb_data_in);
  N_dir_P    = max(dir_P_size);

  for idx_P = 1 : N_dir_P
    % Read data
    % theta_new_read = theta_new_data(idx_P, idx_A);
    % theta_perturb_read = theta_perturb_data(idx_P, idx_A);

    theta_perturb_plot = [theta_perturb_plot, NaN, theta_perturb_data_in{idx_P, idx_A_in}];
    theta_new_plot     = [theta_new_plot, NaN, theta_new_data_in{idx_P, idx_A_in}];
  end

  % idx1 = 1;
  % idx2 = 2;
  % 
  % theta_perturb_plot = theta_perturb_data{idx1, idx2};
  % theta_new_plot     = theta_new_data{idx1, idx2};
  % theta_old_plot     = theta_old_data(idx1, idx2);
  % A_perturb_plot     = A_perturb_data(idx1, idx2);

  x_plot = theta_perturb_plot;
  y_plot = theta_new_plot;

  % y_mod = mod(theta_new_plot, 1);
  % 
  % dx = abs(diff(y_mod));
  % threshold = 0.5;  % You can adjust this if needed
  % breaks = find(dx > threshold);
  % 
  % % Insert NaNs to break the line at discontinuities
  % y_plot = y_mod;
  % x_plot = theta_perturb_plot;
  % for i = length(breaks):-1:1
  %     idx = breaks(i) + 1;
  %     y_plot = [y_plot(1:idx-1), NaN, y_plot(idx:end)];
  %     x_plot = [x_plot(1:idx-1), NaN, x_plot(idx:end)];
  % end

  x_out = x_plot;
  y_out = y_plot;
end