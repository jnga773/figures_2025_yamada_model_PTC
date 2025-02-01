close all; clear all; clc

%%
%----------------------------%
%     Read Data: Folders     %
%----------------------------%
% Run name string identifier
run_name = 'run09_phase_reset_PTC_scan';
% Folder name
dir_data = sprintf('./data/%s/', run_name);
% List all directories
dirs = dir(dir_data);
% Remove ./ and ../
dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dir_sub = {dirs.name};

%------------------------------------%
%     Cycle Through Data Folders     %
%------------------------------------%
% Empty array for A_perturb values
A_perturb = zeros(length(dir_sub), 1);
theta_old = cell(1, length(dir_sub));
theta_new = cell(1, length(dir_sub));

for i = 1 : length(dir_sub)
% for i = 1 : 1
  % Run name
  sub_run_name = {run_name, dir_sub{i}};

  % Bifurcation data
  bd_read = coco_bd_read(sub_run_name);

  % Read A_perturb value
  A_perturb_read = coco_bd_val(bd_read, 1, 'A_perturb');
  % Update array
  A_perturb(i) = A_perturb_read;

  % Read PTC data
  theta_old_read = coco_bd_col(bd_read, 'theta_old');
  theta_new_read = coco_bd_col(bd_read, 'theta_new');

  % Get unique values
  [~, unique_idx] = unique(theta_old_read);
  theta_old_unique = theta_old_read(unique_idx);
  theta_new_unique = theta_new_read(unique_idx);

  % Split theta_old_read into theta_old < 1 and theta_old > 1
  theta_old_lt1 = theta_old_unique(theta_old_unique <= 1.0);
  theta_old_gt1 = theta_old_unique(theta_old_unique > 1.0) - 1.0;

  % Split theta_new_read into theta_old < 1 and theta_old > 1
  theta_new_lt1 = theta_new_unique(theta_old_unique <= 1.0);
  theta_new_gt1 = theta_new_unique(theta_old_unique > 1.0);

  % Check if theta_old_lt1 and theta_old_gt1 cover entire 0 -> 1 range
  lt1_check = false;
  gt1_check = false;

  if round(min(theta_old_lt1), 3) == 0.0 && round(max(theta_old_lt1), 3) == 1.0
    lt1_check = true;
  end
  if round(min(theta_old_gt1), 3) == 0.0 && round(max(theta_old_gt1), 3) == 1.0
    gt1_check = true;
  end

  % If both are true, save the lt1 data
  if lt1_check && gt1_check
    theta_old{i} = theta_old_lt1;
    theta_new{i} = theta_new_lt1;
  end

  % If only one is true, save that data
  if lt1_check && ~gt1_check
    theta_old{i} = theta_old_lt1;
    theta_new{i} = theta_new_lt1;
  elseif ~lt1_check && gt1_check
    theta_old{i} = theta_old_gt1;
    theta_new{i} = theta_new_gt1;
  end

  % If both are false, merge the two together with a NaN in between
  if ~lt1_check && ~gt1_check
    % Check if A_perturb < 0.55
    if A_perturb < 0.55
      % Move theta_new down by one
      theta_new_gt1 = theta_new_gt1 - 1.0;
    end
    
    % Append data
    theta_old{i} = [theta_old_gt1, NaN, theta_old_lt1];
    theta_new{i} = [theta_new_gt1, NaN, theta_new_lt1];
  end
end

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
ax.FontSize = 14;

%--------------------%
%     Plot: PTCs     %
%--------------------%
% Linewidth
lw = 1.5;

% Plot
hold(ax, 'on');

% Plotting index
plot_idx = 7;
plot(ax, theta_old{plot_idx}, theta_new{plot_idx}, LineStyle='-');

% for plot_idx = 1 : length(A_perturb)
%   plot(ax, theta_old{plot_idx}, theta_new{plot_idx}, LineStyle='-');
% end


% Vertical line at intersection point
xline(ax, 0.3176, Color=[0.0 0.0 0.0], LineStyle='-', LineWidth=1.0);

% Plot diagonal line% Plot diagonal line
plot(ax, [0, 1], [0, 1], LineStyle='-', Color=colours(3, :), LineWidth=1.5, ...
     HandleVisibility='off');

% Shade fundamental domain
patch([0, 1, 1, 0], [0, 0, 1, 1], colours(3, :), FaceAlpha=0.2, ...
      EdgeColor='none', HandleVisibility='off')

hold(ax, 'off');

%---------------------%
%     Axis Labels     %
%---------------------%
xlabel(ax, '$\theta_{\mathrm{o}}$');
ylabel(ax, '$\theta_{\mathrm{n}}$');
title_str = sprintf('$A_{\\mathrm{p}} = %.4f$', A_perturb(plot_idx));
title(ax, title_str);

%---------------------%
%     Axis Limits     %
%---------------------%
ax.XAxis.Limits = [0, 1.0];
ax.YAxis.Limits = [-0.25, 2.5];

%----------------------%
%     Figure Stuff     %
%----------------------%
box(ax, 'on');
% grid(ax, 'on');
