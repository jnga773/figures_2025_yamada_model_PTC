clear all; close all;

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
idx_A = 3;

% Cycle through theta_perturb directories and append data
theta_perturb_plot = [];
theta_new_plot     = [];
theta_old_plot     = theta_old_data(1, idx_A);
A_perturb_plot     = A_perturb_data(1, idx_A);

for idx_P = 1 : length(dirs_P)
  % Read data
  % theta_new_read = theta_new_data(idx_P, idx_A);
  % theta_perturb_read = theta_perturb_data(idx_P, idx_A);

  theta_perturb_plot = [theta_perturb_plot, NaN, theta_perturb_data{idx_P, idx_A}];
  theta_new_plot     = [theta_new_plot, NaN, theta_new_data{idx_P, idx_A}];
end

% idx1 = 1;
% idx2 = 2;
% 
% theta_perturb_plot = theta_perturb_data{idx1, idx2};
% theta_new_plot     = theta_new_data{idx1, idx2};
% theta_old_plot     = theta_old_data(idx1, idx2);
% A_perturb_plot     = A_perturb_data(idx1, idx2);

%% TEST PLOT

% x_plot = theta_perturb_plot;
% y_plot = theta_new_plot;

y_mod = mod(theta_new_plot, 1);

dx = abs(diff(y_mod));
threshold = 0.5;  % You can adjust this if needed
breaks = find(dx > threshold);

% Insert NaNs to break the line at discontinuities
y_plot = y_mod;
x_plot = theta_perturb_plot;
for i = length(breaks):-1:1
    idx = breaks(i) + 1;
    y_plot = [y_plot(1:idx-1), NaN, y_plot(idx:end)];
    x_plot = [x_plot(1:idx-1), NaN, x_plot(idx:end)];
end



fprintf('A_perturb = %.4f\n', A_perturb_plot);
fprintf('theta_old = %.4f\n', theta_old_plot);

colours = colororder();

fig = figure(1);
clf;
ax = gca();

width = 6.0;
height = 8.0;

set_figure_dimensions(width, height, scale=2);
% fig.Position = [5, 5, 8, 8];

% Plot
hold(ax, 'on');

% Fundamental domain
patch([-3, 3, 3, -3], [0, 0, 1, 1], colours(3, :), ...
    FaceAlpha=0.2, EdgeColor='none', ...
    HandleVisibility='off');

hold(ax, 'on');

% DTC
% plot(ax, theta_perturb_plot, theta_new_plot, LineStyle='-', Color=colours(1, :));
plot(ax, x_plot, y_plot, LineStyle='-', Color=colours(1, :));
plot(ax, x_plot, y_plot+1, LineStyle='-', Color=colours(1, :));
plot(ax, x_plot, y_plot-1, LineStyle='-', Color=colours(1, :));

plot(ax, x_plot+1, y_plot, LineStyle='-', Color=colours(1, :));
plot(ax, x_plot+1, y_plot+1, LineStyle='-', Color=colours(1, :));
plot(ax, x_plot+1, y_plot-1, LineStyle='-', Color=colours(1, :));

hold(ax, 'off');

% daspect(ax, [1, 1, 1]);
% legend(Interpreter='latex', Location='southwest');

% Make a square
daspect([1, 1, 1]);

% Limits
xlim(ax, [-0.0, 1.0]);
ylim(ax, [-0.5, 1.5]);

% Labels
xlabel(ax, '$\varphi_{\mathrm{d}} / 2 \pi$');
ylabel(ax, '$\vartheta_{\mathrm{n}}$');
% title(ax, title_str);

% Figure stuff
box(ax, 'on');
