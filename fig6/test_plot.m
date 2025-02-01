close all; clear all; clc

%%
%----------------------------%
%     Read Data: Folders     %
%----------------------------%
% Run name string identifier
run_name = 'run10_phase_reset_PTC_scan';
% Folder name
dir_data = sprintf('./data/%s/', run_name);
% List all directories
dirs = dir(dir_data);
% Remove ./ and ../
dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dir_sub = {dirs.name};

%------------------------%
%     Read Data: PTC     %
%------------------------%
% Pick an index to read and plot
read_idx = 7;

% Run name
sub_run_name = {run_name, dir_sub{read_idx}};

% Bifurcation data
bd_read = coco_bd_read(sub_run_name);

% Read A_perturb value
A_perturb_read = coco_bd_val(bd_read, 1, 'A_perturb');

% Read PTC data
theta_old_read = coco_bd_col(bd_read, 'theta_old');
theta_new_read = coco_bd_col(bd_read, 'theta_new');

% Split into theta_old < 1 and theta_old > 1
theta_old_lt1 = theta_old_read(theta_old_read <= 1.0);
theta_old_gt1 = theta_old_read(theta_old_read > 1.0);

theta_new_lt1 = theta_new_read(theta_old_read <= 1.0);
theta_new_gt1 = theta_new_read(theta_old_read > 1.0);

%%
%-------------------%
%     Plot Data     %
%-------------------%
% Default colour order (matplotlib)
colours = colororder();

fig = figure(1); clf;
fig.Name = 'PTC Scan';

% Figure dimensions
fig.Units = 'inches';
fig.Position = [5, 5, 6, 4];

% Figure pdf settings
fig.PaperUnits = fig.Units;
fig.PaperPosition = fig.Position;
fig.PaperSize = fig.Position(3:4);

% Axis setup
tiles = tiledlayout(1, 1, Padding='compact', TileSpacing='compact');
ax = nexttile;
ax.FontSize = 8;

%--------------------%
%     Plot: PTCs     %
%--------------------%
% Linewidth
lw = 1.5;

% Plot
hold(ax, 'on');
plot(ax, theta_old_lt1, theta_new_lt1, LineStyle='-');
plot(ax, theta_old_gt1, theta_new_gt1, LineStyle='--');
hold(ax, 'off');

%-------------------%
%     Hold Axis     %
%-------------------%
hold(ax, 'off');

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 2.0];
ax.YAxis.Limits = [-0.25, 2.0];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');