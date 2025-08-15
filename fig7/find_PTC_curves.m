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

    % theta_old data
    theta_old_out{1, idx} = theta_old_lt1{data_idx};
    theta_old_out{2, idx} = theta_old_gt1{data_idx}-1;
    % theta_new data
    theta_new_out{1, idx} = theta_new_lt1{data_idx};
    theta_new_out{2, idx} = theta_new_gt1{data_idx};
    % A_perturb data
    A_perturb_out{1, idx} = A_perturb(data_idx) * ones(1, length(theta_old_lt1{data_idx}));
    A_perturb_out{2, idx} = A_perturb(data_idx) * ones(1, length(theta_old_gt1{data_idx}));

  end
end