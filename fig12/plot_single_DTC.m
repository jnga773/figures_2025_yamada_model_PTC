%-------------------%
%     Read Data     %
%-------------------%
run_base = 'run07_phase_reset_DTC_scan';

run_sub = {'run_01', 'run_02', 'run_03'};

% Cycle through and read data
theta_old     = zeros(1, length(run_sub));
theta_new     = cell(1, length(run_sub));
A_perturb     = zeros(1, length(run_sub));
theta_perturb = cell(1, length(run_sub));

for idx = 1 : length(run_sub)
  run_in = {run_base, run_sub{idx}};
  bd_read = coco_bd_read(run_in);

  % theta_old
  theta_old_read     = coco_bd_col(bd_read, 'theta_old');
  % theta_new
  theta_new_read     = coco_bd_col(bd_read, 'theta_new');
  % A_perturb
  A_perturb_read     = coco_bd_col(bd_read, 'A_perturb');
  % theta_perturb
  theta_perturb_read = coco_bd_col(bd_read, 'theta_perturb');

  % Save data
  theta_old(idx)     = theta_old_read(1);
  theta_new{idx}     = theta_new_read;
  A_perturb(idx)     = A_perturb_read(1);
  theta_perturb{idx} = theta_perturb_read;
end

%%
%-------------------%
%     Plot Data     %
%-------------------%
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
% plot(ax, theta_perturb, theta_new, LineStyle='-', Color=colours(1, :));

for idx = 1 : length(run_sub)
  theta_perturb_plot = theta_perturb{idx};
  theta_new_plot     = theta_new{idx};

  label_str = sprintf('$A = %.4f$', A_perturb(idx));

  plot(ax, theta_perturb_plot, theta_new_plot, ...
       LineStyle='-', Color=colours(idx, :), ...
       DisplayName=label_str);
end

% daspect(ax, [1, 1, 1]);
legend(Interpreter='latex', Location='southwest');

% Limits
xlim(ax, [0.0, 0.5]);
ylim(ax, [-1.5, 2.0]);

% Labels
xlabel(ax, '$\varphi_{\mathrm{d}} / 2 \pi$');
ylabel(ax, '$\vartheta_{\mathrm{n}}$');
% title_str = sprintf('$\\vartheta_{\\mathrm{o}} = %.4f, A_{\\mathrm{p}} = %.4f$', theta_old(1), A_perturb(1));
title_str = sprintf('$\\vartheta_{\\mathrm{o}} = %.4f$', theta_old(1));
title(ax, title_str);

% Figure stuff
box(ax, 'on');
