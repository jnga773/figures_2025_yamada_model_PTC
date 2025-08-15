function [theta_old_out, A_perturb_out, theta_new_out] = pad_data(data_in, MX_check, N_mesh)
  % [theta_old_out, A_perturb_out, theta_new_out] = pad_data(data_in)
  %
  % Read and pad the "before hole" data

  arguments
    data_in struct
    MX_check = '';
    N_mesh double = 1000;
  end

  % Read data
  theta_old_in = data_in.theta_old;
  theta_new_in = data_in.theta_new;
  A_perturb_in = data_in.A_perturb;

  % Empty arrays for max data length, and max and min theta_old
  array_theta_new_start = zeros(1, length(theta_old_in));
  array_theta_new_end   = zeros(1, length(theta_old_in));
  
  % Cycle through runs
  for idx = 1 : length(theta_old_in)
    % Sort data
    [~, theta_new_sort] = sort_data(theta_old_in{idx}, theta_new_in{idx});

    % Get theta_new end points
    array_theta_new_start(idx) = theta_new_sort(1);
    array_theta_new_end(idx)   = theta_new_sort(end);

  end

  % Find min and max theta_new_values
  if strcmp(MX_check, 'gt1')
    [theta_new_start, ~] = min(array_theta_new_start);
    [theta_new_end, ~]   = max(array_theta_new_end);
  elseif strcmp(MX_check, 'lt1')
    [theta_new_start, ~] = max(array_theta_new_start);
    [theta_new_end, ~]   = min(array_theta_new_end);
  end

  % Output data
  theta_old_out = [];
  theta_new_out = [];
  A_perturb_out = [];

  % Cycle through all data arrays and interpolate data
  for idx = 1 : length(theta_old_in)
    % Read temp data
    theta_old_read = theta_old_in{idx};
    theta_new_read = theta_new_in{idx};
    A_perturb_read = A_perturb_in(idx);

    if isempty(MX_check)
      theta_old_interp = linspace(0.0, 1.0, N_mesh);
    elseif strcmp(MX_check, 'lt1')
      theta_old_interp = linspace(theta_old_read(1), 1.0, round(0.5 * N_mesh));
    elseif strcmp(MX_check, 'gt1')
      theta_old_interp = linspace(0.0, theta_old_read(end), round(0.5 * N_mesh));
    end

    % Get unique indices
    [~, unique_idx] = unique(theta_old_read);
    theta_old_read  = theta_old_read(unique_idx);
    theta_new_read  = theta_new_read(unique_idx);
      
    % Interpolate data
    theta_new_interp = interp1(theta_old_read, theta_new_read, theta_old_interp);
    A_perturb_interp = A_perturb_read * ones(1, length(theta_old_interp));

    % Change max value for holes
    % if idx ~= 1 || idx ~= length(theta_old_in)
       if strcmp(MX_check, 'lt1')
        % Set first value to max
        % theta_new_interp(1)   = theta_new_start;
        theta_new_interp(1)   = array_theta_new_start(idx);
       elseif strcmp(MX_check, 'gt1')
        % theta_new_interp(end) = theta_new_end;
        theta_new_interp(end) = array_theta_new_end(idx);
      end
    % end

    % Append the data
    theta_old_out = [theta_old_out; theta_old_interp];
    theta_new_out = [theta_new_out; theta_new_interp];
    A_perturb_out = [A_perturb_out; A_perturb_interp];

  end
end

