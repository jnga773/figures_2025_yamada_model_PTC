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
