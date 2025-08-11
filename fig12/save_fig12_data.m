function save_fig12_data(run_in, filename_in)
  % save_fig12_data(run_in, filename_in)
  %
  % Reads DTC scan data from run_in and saves data to filename_in
  %
  % Parameters
  % ----------
  % run_in : char
  %    String identifier for the COCO run data.
  % filename_in : char
  %    Name of output data file to save data to.
  %
  % See Also
  % --------
  % coco_bd_read, coco_bd_col, coll_read_solution

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    filename_in char
  end

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
  A_perturb_data     = zeros(length(dirs_P), 3);
  theta_old_data     = zeros(length(dirs_P), 3);
  theta_new_data     = cell(length(dirs_P), 3);
  theta_perturb_data = cell(length(dirs_P), 3);

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
      theta_new_data{idx_P, idx_DTC}     = theta_new_read+1;
      theta_perturb_data{idx_P, idx_DTC} = theta_perturb_read;
    end
  end

  %-------------------%
  %     Sort Data     %
  %-------------------%
  % Empty arrays for output data
  A_perturb     = zeros(3, 1);
  theta_old     = zeros(3, 1);
  theta_new     = cell(3, 1);
  theta_perturb = cell(3, 1);

  % Cycle through A_perturb data runs and sort data
  for idx_A = 1 : 3
    % Parameters
    A_perturb(idx_A) = max(A_perturb_data(:, idx_A));
    theta_old(idx_A) = max(theta_old_data(:, idx_A));

    % Sort data
    [theta_p_temp, theta_n_temp] = combine_data_runs(theta_perturb_data, theta_new_data, idx_A);

    % DTC data
    theta_new{idx_A}     = theta_n_temp;
    theta_perturb{idx_A} = theta_p_temp;
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

function [x_out, y_out] = combine_data_runs(theta_perturb_data_in, theta_new_data_in, idx_A_in)
  % [x_out, y_out] = combine_data_runs(theta_perturb_data_in, theta_new_data_in, idx_A_in)
  %
  % Cycle through theta_perturb directories and append data and also maybe mod-1
  % it maybe if you want :)
  %
  % Parameters
  % ----------
  % theta_perturb_data_in : struct
  %    Structure containing the theta_perturb scan data.
  % theta_new_data_in : struct
  %    Structure containing the theta_new scan data.
  % idx_A : double
  %    Index for a specific perturbation amplitude run.
  %
  % Returns
  % -------
  % x_out, y_out : struct
  %    Data which has been sorted, and modulo-ed 1-ed.
  %
  % See Also
  % --------
  % size, max, mod

  theta_perturb_plot = [];
  theta_new_plot     = [];

  %--------------------------------------------------%
  %     Combine Different theta_perturb Run Data     %
  %--------------------------------------------------%
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

  x_plot = theta_perturb_plot;
  y_plot = theta_new_plot;

  %----------------------------%
  %     Take Data Modulo-1     %
  %----------------------------%
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

  %----------------%
  %     Output     %
  %----------------%
  x_out = x_plot;
  y_out = y_plot;

end