function [theta_old_out, theta_perturb_out, A_perturb_out] = calc_interpolated_data()
  % [theta_old_out, theta_perturb_out, A_perturb_out] = calc_interpolated_data()
  %
  % Interpolates theta_old, theta_perturb, and A_perturb over a finer
  % range.
  %-------------------------------------------------------------------------%
  %                              Calculate Data                             %
  %-------------------------------------------------------------------------%
  % Calculate data
  [theta_old, theta_perturb, A_perturb] = calc_intersection_points();

  % Normalise theta_perturb by pi
  theta_perturb = theta_perturb / pi;

  %-----------------------------------------%
  %     Sort Data: Increasing theta_old     %
  %-----------------------------------------%
  % Sort
  [~, sort_idx] = sort(theta_old);

  % Sort by theta_perturb
  theta_old     = theta_old(sort_idx);
  theta_perturb = theta_perturb(sort_idx);
  A_perturb     = A_perturb(sort_idx);

  %---------------------------------------------%
  %     Sort Data: Increasing theta_perturb     %
  %---------------------------------------------%
  % % Sort
  % [~, sort_idx] = sort(theta_perturb);
  % 
  % % Sort by theta_perturb
  % theta_old     = theta_old(sort_idx);
  % theta_perturb = theta_perturb(sort_idx);
  % A_perturb     = A_perturb(sort_idx);

  %------------------------------------------%
  %     Shift Data: By max theta_perturb     %
  %------------------------------------------%
  % Find max point of theta_old
  [max_val, max_idx] = max(theta_perturb);

  % Shift data around
  theta_perturb = [theta_perturb(1:max_idx)-2; theta_perturb(max_idx+1:end)];
  theta_old     = [theta_old(1:max_idx); theta_old(max_idx+1:end)];
  A_perturb     = [A_perturb(1:max_idx); A_perturb(max_idx+1:end)];

  %---------------------%
  %     Extend Data     %
  %---------------------%
  theta_perturb = [theta_perturb-2; theta_perturb; theta_perturb+2];
  theta_old     = [theta_old-1; theta_old; theta_old+1];
  A_perturb     = [A_perturb; A_perturb; A_perturb];

  %%
  %--------------------------------------------%
  %     Interpolate Data: By theta_perturb     %
  %--------------------------------------------%
  % Ranges
  range_min = -1.0;
  range_max = 2.5;
  drange    = 1e-3;

  % Get data over range theta_perturb in [-0.5, 1.0]
  range_idx = theta_perturb >= range_min & theta_perturb <= range_max;

  % Make data smaller range
  theta_perturb_small_range = theta_perturb(range_idx);
  theta_old_small_range     = theta_old(range_idx);
  A_perturb_small_range     = A_perturb(range_idx);

  [~, unique_idx] = unique(theta_perturb_small_range);
  theta_perturb_small_range = theta_perturb_small_range(unique_idx);
  theta_old_small_range = theta_old_small_range(unique_idx);
  A_perturb_small_range = A_perturb_small_range(unique_idx);

  % Smooth theta_perturb data to interpolate over
  theta_perturb_interpolate = range_min : drange : range_max;
  theta_perturb_interpolate = theta_perturb_interpolate';
  % theta_perturb_interpolate = unique(theta_perturb_interpolate);

  % Inteprolate data
  theta_old_interpolate = interp1(theta_perturb_small_range, theta_old_small_range, theta_perturb_interpolate);
  A_perturb_interpolate = interp1(theta_perturb_small_range, A_perturb_small_range, theta_perturb_interpolate);

  % Halve theta_perturb
  theta_perturb_interpolate = 0.5 * theta_perturb_interpolate;

  %-----------------------------%
  %     Print Intersections     %
  %-----------------------------%
  % Print values at theta_perturb = 0 and "pi"
  idx_G    = find(theta_perturb_interpolate == 0.0);
  idx_I    = find(theta_perturb_interpolate == 0.25);
  % idx_half = find(theta_perturb_interpolate == 0.125);

  perturb_angle = 45;
  idx_half = find(round(theta_perturb_interpolate, 3) == round(deg2rad(perturb_angle) / (2 * pi), 3));
  idx_half = idx_half(1);

  % Intersection for G perturbation
  fprintf('G intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
          mod(theta_old_interpolate(idx_G), 1), A_perturb_interpolate(idx_G));

  % Intersection for I perturbation
  fprintf('I intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
          mod(theta_old_interpolate(idx_I), 1), A_perturb_interpolate(idx_I));

  % Intersection for 45deg perturbation
  fprintf('45deg intersection at theta_old = %.5f at A_perturb = %.5f\n', ...
          mod(theta_old_interpolate(idx_half), 1), A_perturb_interpolate(idx_half));

  %----------------%
  %     Output     %
  %----------------%
  theta_old_out     = theta_old_interpolate;
  theta_perturb_out = theta_perturb_interpolate;
  A_perturb_out     = A_perturb_interpolate;

end