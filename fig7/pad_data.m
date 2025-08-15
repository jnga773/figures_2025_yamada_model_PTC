function [theta_old_out, A_perturb_out, theta_new_out] = pad_data(data_in, MX_check)
  % [theta_old_out, A_perturb_out, theta_new_out] = pad_data(data_in)
  %
  % Read and pad the "before hole" data

  arguments
    data_in
    MX_check = '';
  end

  % Read data
  theta_old_read = data_in.theta_old;
  theta_new_read = data_in.theta_new;
  A_perturb_read = data_in.A_perturb;

  % Empty arrays for max data length, and max and min theta_old
  data_lengths          = zeros(1, length(theta_old_read));
  array_theta_old_start = zeros(1, length(theta_old_read));
  array_theta_old_end   = zeros(1, length(theta_old_read));
  array_theta_new_start = zeros(1, length(theta_old_read));
  array_theta_new_end   = zeros(1, length(theta_old_read));
  
  % Cycle through runs
  for idx = 1 : length(theta_old_read)
    % Sort data
    [theta_old_sort, theta_new_sort] = sort_data(theta_old_read{idx}, theta_new_read{idx});

    % Get data lengths
    data_lengths(idx)          = length(theta_old_sort);
    % Get theta_old end points
    array_theta_old_start(idx) = theta_old_sort(1);
    array_theta_old_end(idx)   = theta_old_sort(end);
    % Get theta_new end points
    array_theta_new_start(idx) = theta_new_sort(1);
    array_theta_new_end(idx)   = theta_new_sort(end);

  end

  % Find min and max theta_old values
  [theta_old_start, ~] = min(array_theta_old_start);
  [theta_old_end, ~]   = max(array_theta_old_end);

  % Find min and max theta_new_values
  if strcmp(MX_check, 'gt1')
    [theta_new_start, ~] = min(array_theta_new_start);
    [theta_new_end, ~]   = max(array_theta_new_end);
  elseif strcmp(MX_check, 'lt1')
    [theta_new_start, ~] = max(array_theta_new_start);
    [theta_new_end, ~]   = min(array_theta_new_end);
  end

  % Find the data array with the maximum length
  [~, max_idx] = max(data_lengths);
  
  % Get the maximum length data arrays
  theta_old_max_length = theta_old_read{max_idx};

  % Append min amd max values either side
  theta_old_max_length = [theta_old_start, theta_old_max_length, theta_old_end];

  % Output data
  theta_old_out = [];
  theta_new_out = [];
  A_perturb_out = [];

  % Cycle through all data arrays and interpolate data
  for idx = 1 : length(theta_new_read)
    % Read temp data
    theta_old_temp = theta_old_read{idx};
    theta_new_temp = theta_new_read{idx};
    A_perturb_temp = A_perturb_read(idx);

    if isempty(MX_check)
      % Append min and max values either side
      theta_old_temp = [theta_old_start, theta_old_temp, theta_old_end];
      theta_new_temp = [theta_new_temp(1), theta_new_temp, theta_new_temp(end)];

    elseif strcmp(MX_check, 'lt1')
      % Append max value at the end
      theta_old_temp = [theta_old_temp, theta_old_end];
      theta_new_temp = [theta_new_temp, array_theta_new_end(idx)];

    elseif strcmp(MX_check, 'gt1')
      % Append min value at the start
      theta_old_temp = [theta_old_start, theta_old_temp];
      theta_new_temp = [array_theta_new_start(idx), theta_new_temp];
    end

    % Get unique indices
    [~, unique_idx] = unique(theta_old_temp);
    theta_old_temp  = theta_old_temp(unique_idx);
    theta_new_temp  = theta_new_temp(unique_idx);
      
    % Interpolate data
    theta_new_interp = interp1(theta_old_temp, theta_new_temp, theta_old_max_length);
    % The "interpolated" theta_old array is just the max length array, and the
    % A_perturb array is just a ones array with the same length.
    theta_old_interp = theta_old_max_length;
    A_perturb_interp = A_perturb_temp * ones(1, length(theta_old_interp));

    % Check for NaNs
    if strcmp(MX_check, 'lt1')
      % Find NaNs
      nan_idx = isnan(theta_new_interp);

      % If 'lt1' data, check for NaNs in theta_old < theta_old_min(idx)
      lt1_idx = theta_old_max_length < array_theta_old_start(idx);

      % Change values of theta_old_interp to min value
      theta_old_interp(nan_idx(lt1_idx)) = array_theta_old_start(idx);

      % Set NaNs of lt1_idx to min val
      theta_new_interp(nan_idx(lt1_idx)) = theta_new_start;

    elseif strcmp(MX_check, 'gt1')
      % If 'gt1' data, check for NaNs in theta_old > theta_old_max(idx)
      gt1_idx = theta_old_max_length > array_theta_old_end(idx);

      % Change values of theta_old_interp to min value
      theta_old_interp(gt1_idx) = array_theta_old_end(idx);

      % Set NaNs of lt1_idx to min val
      theta_new_interp(gt1_idx) = theta_new_end;

      % % Find nans
      % nan_idx = isnan(theta_new_interp);
      % theta_new_interp(nan_idx) = array_theta_new_start(idx);

    end

    % Append the data
    theta_old_out = [theta_old_out; theta_old_interp];
    theta_new_out = [theta_new_out; theta_new_interp];
    A_perturb_out = [A_perturb_out; A_perturb_interp];

  end
end

