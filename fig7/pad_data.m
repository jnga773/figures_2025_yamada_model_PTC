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
  data_lengths  = [];
  max_theta_olds = [];
  max_theta_news = [];
  min_theta_olds = [];
  min_theta_news = [];
  
  % Cycle through runs
  for idx = 1 : length(theta_old_read)
    % Sort data
    [theta_old_sort, theta_new_sort] = sort_data(theta_old_read{idx}, theta_new_read{idx});

    % Get data lengths
    data_lengths = [data_lengths, length(theta_old_sort)];

    % Get end values
    min_theta_olds = [min_theta_olds, theta_old_sort(1)];
    min_theta_news = [min_theta_news, theta_new_sort(1)];
    max_theta_olds = [max_theta_olds, theta_old_sort(end)];
    max_theta_news = [max_theta_news, theta_new_sort(end)];

  end

  % Find the max and min theta_old values
  [theta_old_min, theta_old_min_idx] = min(min_theta_olds);
  [theta_old_max, theta_old_max_idx] = max(max_theta_olds);
  [theta_new_min, theta_new_min_idx] = min(min_theta_news);
  [theta_new_max, theta_new_max_idx] = max(min_theta_news);

  % Find the data array with the maximum length
  [~, max_idx] = max(data_lengths);
  
  % Get the maximum length data arrays
  theta_old_max_length = theta_old_read{max_idx};

  % Append min amd max values either side
  theta_old_max_length = [theta_old_min, theta_old_max_length, theta_old_max];

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
      theta_old_temp = [theta_old_min, theta_old_temp, theta_old_max];
      theta_new_temp = [theta_new_temp(1), theta_new_temp, theta_new_temp(end)];

    elseif strcmp(MX_check, 'lt1')
      % Append max value at the end
      theta_old_temp = [theta_old_temp, theta_old_max];
      theta_new_temp = [theta_new_temp, theta_new_temp(end)];

    elseif strcmp(MX_check, 'gt1')
      % Append min value at the start
      theta_old_temp = [theta_old_min, theta_old_temp];
      theta_new_temp = [theta_new_temp(1), theta_new_temp];
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
      lt1_idx = theta_old_max_length < min_theta_olds(idx);

      % Change values of theta_old_interp to min value
      theta_old_interp(nan_idx(lt1_idx)) = min_theta_olds(idx);

      % Set NaNs of lt1_idx to min val
      theta_new_interp(nan_idx(lt1_idx)) = min_theta_news(idx);

    elseif strcmp(MX_check, 'gt1')
      % If 'gt1' data, check for NaNs in theta_old > theta_old_max(idx)
      gt1_idx = theta_old_max_length > max_theta_olds(idx);

      % Change values of theta_old_interp to min value
      theta_old_interp(gt1_idx) = max_theta_olds(idx);

      % Set NaNs of lt1_idx to min val
      theta_new_interp(gt1_idx) = max_theta_news(idx);

      % Find nans
      nan_idx = isnan(theta_new_interp);
      theta_new_interp(nan_idx) = min_theta_news(idx);

    end

    % Append the data
    theta_old_out = [theta_old_out; theta_old_interp];
    theta_new_out = [theta_new_out; theta_new_interp];
    A_perturb_out = [A_perturb_out; A_perturb_interp];

  end
end

function [theta_old_out, theta_new_out] = sort_data(theta_old_in, theta_new_in)
  % 
  %
  % Sorts the data so it's all even and nice, just in cases

  % Get sorted indices
  [~, sort_idx] = sort(theta_old_in);

  % Sort it
  theta_old_out = theta_old_in(sort_idx);
  theta_new_out = theta_new_in(sort_idx);

end
