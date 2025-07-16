clear all; close all;

%% READ DATA
run_in = 'run07_PR_DTC_scan';

%------------------------------%
%     Read Sub-Directories     %
%------------------------------%
% Folder name
dirs_A = sprintf('./data/%s/', run_in);
% List all directories
dirs_A = dir(dirs_A);
% Remove ./ and ../
dirs_A = dirs_A(~ismember({dirs_A.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dirs_A = {dirs_A.name};

%------------------------------------%
%     Cycle Through and Read DTC     %
%------------------------------------%
% Empty cells for data
A_perturb_data     = zeros(length(dirs_A), 2);
theta_old_data     = zeros(length(dirs_A), 2);
theta_new_data     = cell(length(dirs_A), 2);
theta_perturb_data = cell(length(dirs_A), 2);

% Cycle through A_perturb sub-directories
for idx_A = 1 : length(dirs_A)
  % Run name
  A_sub_run = {run_in, dirs_A{idx_A}};

  % Read folders inside dirs_A
  dirs_T = dir(sprintf('./data/%s/%s/', run_in, dirs_A{idx_A}));
  dirs_T = dirs_T(~ismember({dirs_T.name}, {'.', '..', '.DS_Store'}));
  dirs_T = {dirs_T.name};

  % Cycle through theta directories
  for idx_T = 1 : length(dirs_T)
    % Run name
    T_sub_run = {run_in, dirs_A{idx_A}, dirs_T{idx_T}};

    % Read bifurcation data
    bd_read = coco_bd_read(T_sub_run);

    % Get values
    A_perturb_read     = coco_bd_val(bd_read, 1, 'A_perturb');
    theta_old_read     = coco_bd_val(bd_read, 1, 'theta_old');
    theta_new_read     = coco_bd_col(bd_read, 'theta_new');
    theta_perturb_read = coco_bd_col(bd_read, 'theta_perturb');
    
    % Update arrays
    A_perturb_data(idx_A, idx_T)     = A_perturb_read;
    theta_old_data(idx_A, idx_T)     = theta_old_read;
    theta_new_data{idx_A, idx_T}    = theta_new_read;
    theta_perturb_data{idx_A, idx_T} = theta_perturb_read;

  end
end

%% MERGE DATA SET 1
% Get data for first A_perturb value
A1  = A_perturb_data(1, 1);
TO1 = theta_old_data(1, 1);
TN1 = [theta_new_data{1, 1}, theta_new_data{1, 2}]';
TP1 = [theta_perturb_data{1, 1}, theta_perturb_data{1, 2}]';

% Sort
[~, sort_idx] = sort(TP1);
TN1 = TN1(sort_idx);
TP1 = TP1(sort_idx);

%% MERGE DATA SET 2
% Get data for first A_perturb value
A2  = A_perturb_data(2, 1);
TO2 = theta_old_data(2, 1);
TN2 = [theta_new_data{2, 1}]';
TP2 = [theta_perturb_data{2, 1}]';

% Sort
[~, sort_idx] = sort(TP2);
TN2 = TN2(sort_idx);
TP2 = TP2(sort_idx);

%% MERGE DATA SET 3
% Get data for first A_perturb value
A3  = A_perturb_data(3, 1);
TO3 = theta_old_data(3, 1);
TN3 = [theta_new_data{3, 1}, theta_new_data{3, 2}]';
TP3 = [theta_perturb_data{3, 1}, theta_perturb_data{3, 2}]';

% Sort
[~, sort_idx] = sort(TP3);
TN3 = TN3(sort_idx);
TP3 = TP3(sort_idx);

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

hold(ax, 'on');

% DTC
% theta_perturb_plot = theta_perturb{1, 2};
% theta_new_plot     = theta_new{1, 2};
theta_perturb_plot = TP1;
theta_new_plot     = TN1;
plot(ax, theta_perturb_plot, theta_new_plot, LineStyle='-', Color=colours(1, :));

hold(ax, 'off');

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
