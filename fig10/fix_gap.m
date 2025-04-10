function data_out = fix_gap(data_before_hole, data_hole_lt1)
  % Fills in the gap between data_hole_lt1 and the level below of
  % data_before_hole.

  % Get data from end of data_before_hole
  A_perturb_before_hole = data_before_hole.A_perturb(end);
  theta_old_before_hole = data_before_hole.theta_old{end};
  theta_new_before_hole = data_before_hole.theta_new{end} - 1.0;

  % Get data from start of data_hole_lt1
  A_perturb_hole = data_hole_lt1.A_perturb(1);
  theta_old_hole = data_hole_lt1.theta_old{1};
  theta_new_hole = data_hole_lt1.theta_new{1};

  % Find nearest point in theta_old_before_hole to the start point of
  % theta_old_hole
  [min_val, min_idx] = min(abs(theta_old_before_hole - theta_old_hole(1)));

  % Append the data
  A_perturb_out = [A_perturb_before_hole, data_hole_lt1.A_perturb(2:end)];
  theta_old_out = {theta_old_before_hole(min_idx:end), data_hole_lt1.theta_old{2:end}};
  theta_new_out = {theta_new_before_hole(min_idx:end), data_hole_lt1.theta_new{2:end}};

  data_out.A_perturb = A_perturb_out;
  data_out.theta_old = theta_old_out;
  data_out.theta_new = theta_new_out;


end