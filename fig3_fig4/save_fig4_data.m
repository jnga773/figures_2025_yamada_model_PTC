function save_fig4_data(run_in, filename_in)
  % save_fig4_data(run_in, filename_in)
  %
  % Reads PTC data from run_in and saves data to 'filename_in'.

  %--------------------------------%
  %     Read Data: Phase Reset     %
  %--------------------------------%
  % Bifurcation data
  bd_PR = coco_bd_read(run_in);

  % Get solution labels
  label_PR = coco_bd_labs(bd_PR, 'SP');

  % Get theta_old values
  theta_old_run1 = coco_bd_val(bd_PR, label_PR(1), 'theta_old');
  theta_old_run2 = coco_bd_val(bd_PR, label_PR(2), 'theta_old');

  % Get A_perturb value
  A_perturb = coco_bd_val(bd_PR, label_PR(1), 'A_perturb');

  % Read theta_old and theta_new values from bifurcation data
  theta_old_plot = coco_bd_col(bd_PR, 'theta_old');
  theta_new_plot = coco_bd_col(bd_PR, 'theta_new');

  %-------------------%
  %     Save Data     %
  %-------------------%
  % Parameters
  data_out.theta_old_run1 = theta_old_run1;
  data_out.theta_old_run2 = theta_old_run2;
  data_out.A_perturb      = A_perturb;

  % Plotting data
  data_out.theta_old_plot = theta_old_plot;
  data_out.theta_new_plot = theta_new_plot;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end