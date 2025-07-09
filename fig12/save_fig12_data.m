function save_fig8_data(run_in, filename_in)
  % save_fig8_data(run_in, filename_in)
  %
  % Reads DTC scan data from run_in and saves data to filename_in

  %------------------------%
  %     Read Data: DTCs    %
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
  theta_old     = zeros(1, length(dir_sub));
  A_perturb     = zeros(1, length(dir_sub));
  theta_perturb = cell(1, length(dir_sub));
  theta_new     = cell(1, length(dir_sub));

  % Cycle through data sub directories
  for idx = 1 : length(dir_sub)
    % Run name
    sub_run_name = {run_in, dir_sub{idx}};

    % Bifurcation data
    bd_read = coco_bd_read(sub_run_name);

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
    theta_new{idx}     = theta_new_read - 1.0;
    A_perturb(idx)     = A_perturb_read(1);
    theta_perturb{idx} = theta_perturb_read;
  end

  %-------------------%
  %     Save Data     %
  %-------------------%
  % Parameters
  data_out.A_perturb     = A_perturb;
  data_out.theta_old     = theta_old;

  % Plotting data
  data_out.theta_new     = theta_new;
  data_out.theta_perturb = theta_perturb;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end