function save_fig6_data(run_in, filename_in)
  % save_fig6_data(run_in, filename_in)
  %
  % Reads PTC scan data from run_in and saves data to filename_in

  %-----------------------------------%
  %     Read Data: Periodic Orbit     %
  %-----------------------------------%

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

  % Get list of data
  dir_sub_plot = dir_sub;

  % Empty arrays for theta_old and theta_new
  theta_old_data = cell([length(dir_sub_plot), 1]);
  theta_new_data = cell([length(dir_sub_plot), 1]);
  A_perturb_data = zeros([length(dir_sub_plot), 1]);

  % Cycle through each data directory
  for i = 1 : length(dir_sub_plot)
    % Sub folder name
    dir_read = {PTC_run_in, dir_sub_plot{i}};
    fprintf('run_dir:  {%s, %s} \n', dir_read{1}, dir_read{2});

    % Bifurcation data
    bd_read = coco_bd_read(dir_read);

    % Read theta_old and theta_new
    theta_old_read = coco_bd_col(bd_read, 'theta_old');
    theta_new_read = coco_bd_col(bd_read, 'theta_new');
    A_perturb_read = coco_bd_col(bd_read, 'A_perturb');

    % Append to arrays
    theta_old_data{i} = theta_old_read;
    theta_new_data{i} = theta_new_read;
    A_perturb_data(i) = A_perturb_read(1);

  end

  % Read parameters
  theta_perturb = coco_bd_val(bd_read, 1, 'theta_perturb');

  %-------------------%
  %     Save Data     %
  %-------------------%
  % Parameters
  data_out.theta_perturb = theta_perturb;

  % PTC data
  data_out.theta_old     = theta_old_data;
  data_out.theta_new     = theta_new_data;
  data_out.A_perturb     = A_perturb;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end