function [theta_old_out, A_perturb_out, theta_new_out] = find_PTC_curves(filename_in, A_perturb_values)
  % [theta_old_out, A_perturb_out, theta_new_out] = find_PTC_curves(filename_in, A_perturb_values)
  %
  % Finds the PTC curves for a given set of perturbation amplitudes from the data
  % in filename_in.

  % Load data
  load(filename_in, 'A_perturb', 'theta_old_lt1', 'theta_new_lt1', ...
       'theta_old_gt1', 'theta_new_gt1');

  % Empty cells for plotting data
  theta_old_out = cell(2, length(A_perturb_values));
  theta_new_out = cell(2, length(A_perturb_values));
  A_perturb_out = cell(2, length(A_perturb_values));

  % Cycle through A_perturb_values, find the appropriate data index, and
  % save data
  for idx = 1 : length(A_perturb_values)
    % Look for data
    data_idx = find(round(A_perturb, 4) == A_perturb_values(idx));

    % interpolate lt1 data
    [X1, Y1, Z1] = interp_data(theta_old_lt1{data_idx}, A_perturb(data_idx), theta_new_lt1{data_idx});
    % interpolate gt1 data
    [X2, Y2, Z2] = interp_data(theta_old_gt1{data_idx}, A_perturb(data_idx), theta_new_gt1{data_idx});

    % theta_old data
    theta_old_out{1, idx} = X1;
    theta_old_out{2, idx} = X2-1;
    % A_perturb data
    A_perturb_out{1, idx} = Y1;
    A_perturb_out{2, idx} = Y2;
    % theta_new data
    theta_new_out{1, idx} = Z1;
    theta_new_out{2, idx} = Z2;

  end
end

function [x_interp, y_interp, z_interp] = interp_data(x_in, y_in, z_in, N_mesh)
  % [x_out, y_out, z_out] = reinterp_data(x_in, y_in, z_in, N_mesh)
  %
  % Interpolate input data for fewer points

  arguments
    x_in double
    y_in double
    z_in double
    N_mesh double = 500;
  end

  % Get unique data
  [~, unique_idx] = unique(x_in);
  x_in = x_in(unique_idx);
  z_in = z_in(unique_idx);

  % Set x_in as linspace array
  x_interp = linspace(min(x_in), max(x_in), N_mesh);

  % y_in interp is just ones
  y_interp = y_in * ones(1, length(x_interp));

  % Interpolate z_in
  z_interp = interp1(x_in, z_in, x_interp);

end