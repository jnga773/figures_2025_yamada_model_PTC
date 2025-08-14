function [theta_old_out, A_perturb_out, theta_new_out] = find_PTC_curves(filename_in, A_perturb_values)
  % [theta_old_out, A_perturb_out, theta_new_out] = find_PTC_curves(filename_in, A_perturb_values)
  %
  % Finds the PTC curves for a given set of perturbation amplitudes from the data
  % in filename_in.

  % Load data
  load(filename_in, 'A_perturb', 'theta_old_lt1', 'theta_new_lt1', ...
       'theta_old_gt1', 'theta_new_gt1');

  % Find plotting indices
  plot_idx = zeros(length(A_perturb_values), 1);
  for i = 1 : length(A_perturb_values)
    plot_idx(i) = find(round(A_perturb, 4) == A_perturb_values(i));
  end

  % Empty cells for plotting data
  theta_old_out = cell(1, length(plot_idx));
  theta_new_out = cell(1, length(plot_idx));
  A_perturb_out = cell(1, length(plot_idx));

  for i = 1 : length(A_perturb_values)
    % Plot data index
    idx = plot_idx(i);

    % Logical checks
    lt1_check = false;
    gt1_check = false;

    if round(theta_old_lt1{idx}(end) - theta_old_lt1{idx}(1), 3) == 1.0
      lt1_check = true;
    end
    if round(theta_old_gt1{idx}(end) - theta_old_gt1{idx}(1), 3) == 1.0
      gt1_check = true;
    end

    % Grab data
    if lt1_check && gt1_check
      % Both checks are true; take the lt1 data
      theta_old_out{i} = theta_old_lt1{idx};
      theta_new_out{i} = theta_new_lt1{idx};
      A_perturb_out{i} = A_perturb(idx) * ones(1, length(theta_old_lt1{idx}));
    elseif lt1_check && ~gt1_check
      % Only the lt1 check is true
      theta_old_out{i} = theta_old_lt1{idx};
      theta_new_out{i} = theta_new_lt1{idx};
      A_perturb_out{i} = A_perturb(idx) * ones(1, length(theta_old_lt1{idx}));
    elseif ~lt1_check && gt1_check
      % Onle the gt1 check is true
      theta_old_out{i} = theta_old_gt1{idx}-1.0;
      theta_new_out{i} = theta_new_gt1{idx};
      A_perturb_out{i} = A_perturb(idx) * ones(1, length(theta_old_gt1{idx}));
    else
      % Neither of the checks are true
      theta_old_out{i} = [theta_old_gt1{idx}-1.0, nan, theta_old_lt1{idx}];
      theta_new_out{i} = [theta_new_gt1{idx}, nan, theta_new_lt1{idx}];
      A_perturb_out{i} = A_perturb(idx) * ones(1, length([theta_old_gt1{idx}, nan, theta_old_lt1{idx}]));
    end
  end
end