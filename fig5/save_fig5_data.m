function save_fig6_data(run_in, filename_in)
  % save_fig6_data(run_in, filename_in)
  %
  % Reads PTC scan data from run_in and saves data to filename_in
  %
  % Parameters
  % ----------
  % run_names_in : char
  %    Run name string identifier.
  % filename_in : char
  %    Name of output data file to save data to.
  %
  % See Also
  % --------
  % coco_bd_read, coco_bd_labs, coll_read_solution, ep_read_solution

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    filename_in char
  end
  
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
    sub_run_name = {run_in, dir_sub{i}};

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
      % Append data
      theta_old{i} = [theta_old_gt1, NaN, theta_old_lt1];
      theta_new{i} = [theta_new_gt1, NaN, theta_new_lt1];
    end

    % Check if theta_new starts in the fundamental domain. If not, shift it
    if theta_old{i}(1) < 0.0
      theta_new{i} = theta_new{i} + 1.0;
    elseif theta_old{i}(1) > 1.0
      theta_new{i} = theta_new{i} - 1.0;
    end
  end

  %-------------------%
  %     Save Data     %
  %-------------------%
  % PTC data
  data_out.theta_old     = theta_old;
  data_out.theta_new     = theta_new;
  data_out.A_perturb     = A_perturb;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end