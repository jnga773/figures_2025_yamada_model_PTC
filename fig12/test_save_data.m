clear all; close all;

%% READ DATA
run_in = 'run07_phase_reset_DTC_scan';

%------------------------%
%     Read Data: PTCs    %
%------------------------%
% Folder name
dir_data = sprintf('./data/%s/', run_in);
% List all directories
dirs = dir(dir_data);
% Remove ./ and ../
dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dir_sub = {dirs.name};

% Empty cells
theta_old     = zeros(length(dir_sub), 1);
A_perturb     = zeros(length(dir_sub), 1);
% theta_perturb = cell(length(dir_sub), 1);
% theta_new     = cell(length(dir_sub), 1);
SP_labs       = cell(length(dir_sub), 1);
theta_perturb = zeros(length(dir_sub), 1);
theta_new     = zeros(length(dir_sub), 1);

% Cycle through all sub directories and read data
for idx = 1 : length(dir_sub)
  % Run name
  sub_run_name = {run_in, dir_sub{idx}};

  % Bifurcation data
  bd_read = coco_bd_read(sub_run_name);

  % Read SP labels
  labs_read = coco_bd_labs(bd_read, 'SP');
  labs_read = sort(labs_read);

  % Save labels
  SP_labs{idx} = labs_read;

  lab_read = labs_read(1);

  % Read data
  theta_old(idx)     = coco_bd_val(bd_read, lab_read, 'theta_old');
  A_perturb(idx)     = coco_bd_val(bd_read, lab_read, 'A_perturb');
  theta_new(idx)     = coco_bd_val(bd_read, lab_read, 'theta_new');
  theta_perturb(idx) = coco_bd_val(bd_read, lab_read, 'theta_perturb');

end

% Sort data by theta_perturb
[~, sort_idx] = sort(theta_perturb);
theta_perturb = theta_perturb(sort_idx);
theta_old     = theta_old(sort_idx);
A_perturb     = A_perturb(sort_idx);
theta_new     = theta_new(sort_idx);

% Renormalise theta_new
theta_new = theta_new - 1.0;

%% PLOT DATA
colours = colororder();

fig = figure(1); clf;
ax = gca();

width = 6.0;
height = 8.0;

set_figure_dimensions(width, height);
% fig.Position = [5, 5, 8, 8];

% Plot
hold(ax, 'on');

% Fundamental domain
patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
    FaceAlpha=0.2, EdgeColor='none', ...
    HandleVisibility='off');

% DTC
plot(ax, theta_perturb, theta_new, LineStyle='-', Color=colours(1, :));

% daspect(ax, [1, 1, 1]);
legend(Interpreter='latex', Location='southwest');

% Limits
xlim(ax, [0.0, 1.0]);
ylim(ax, [-1.5, 2.0]);

% Labels
xlabel(ax, '$\varphi_{\mathrm{d}} / 2 \pi$');
ylabel(ax, '$\vartheta_{\mathrm{n}}$');
% title(ax, title_str);

% Figure stuff
box(ax, 'on');
