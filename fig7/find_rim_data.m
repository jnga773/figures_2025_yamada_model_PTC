function [theta_old_out, A_perturb_out, theta_new_out] = find_rim_data(X_in, Y_in, Z_in, lt1_gt1)
  % [theta_old_out, A_perturb_out, theta_new_out] = find_rim_data(X_in, Y_in, Z_in)
  %
  % Finds and outputs the PTC data around the rim of the "volcano" to plot.

  % Find column with max theta_new value and output that column
  if strcmp(lt1_gt1, 'gt1')
    theta_old_out = X_in(:, end);
    A_perturb_out = Y_in(:, end);
    theta_new_out = Z_in(:, end);
  elseif strcmp(lt1_gt1, 'lt1')
    theta_old_out = X_in(:, 1);
    A_perturb_out = Y_in(:, 1);
    theta_new_out = Z_in(:, 1);
  end

  % Cut off the edges
  theta_old_out = theta_old_out(2:end-1);
  A_perturb_out = A_perturb_out(2:end-1);
  theta_new_out = theta_new_out(2:end-1);

end