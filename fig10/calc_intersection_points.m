function [theta_old_out, theta_perturb_out, A_perturb_out] = calc_intersection_points()
  %
  %
  % Calculate the directional vector to perturb orbit into the manifold in the
  % G-I plane.

  %-------------------%
  %     Read Data     %
  %-------------------%
  load('../data_files/fig2_data.mat', 'xpos', 'Wq_s', 'xbp_PO', 'tbp_PO', 'T_PO');

  % Normalise time data by the period
  tbp_PO = tbp_PO / T_PO;
  % Rename xbp
  % xbp_PO = xbp_PO;

  %-------------------------------------%
  %     Find Distance from Manifold     %
  %-------------------------------------%
  % Break up components of periodic orbit
  G_PO = xbp_PO(:, 1);
  Q_PO = xbp_PO(:, 2);
  I_PO = xbp_PO(:, 3);

  % Break up components of stable manifold of q
  G_W = Wq_s(:, 1);
  Q_W = Wq_s(:, 2);
  I_W = Wq_s(:, 3);

  % Empty array for differences
  G_diff = zeros(length(G_PO), length(G_W));
  Q_diff = zeros(length(Q_PO), length(Q_W));
  I_diff = zeros(length(I_PO), length(I_W));

  % Cycle through periodic orbit and stable manifold to find minimum
  for i = 1 : length(Q_PO)
    for j = 1 : length(Q_W)
      % Calculate difference
      G_diff(i, j) = abs(G_PO(i) - G_W(j));
      Q_diff(i, j) = abs(Q_PO(i) - Q_W(j));
      I_diff(i, j) = abs(I_PO(i) - I_W(j));

    end
  end

  %----------------------------------------%
  %     Find Closest Point to Manifold     %
  %----------------------------------------%
  % Calculate minimum difference point
  G_min_val = zeros(length(Q_PO), 1); G_min_idx = zeros(length(Q_PO), 1);
  Q_min_val = zeros(length(Q_PO), 1); Q_min_idx = zeros(length(Q_PO), 1);
  I_min_val = zeros(length(Q_PO), 1); I_min_idx = zeros(length(Q_PO), 1);

  % Cycle through differences
  for i = 1 : length(Q_PO)
    % Calculate minimum
    [G_min_temp, G_idx_temp] = min(G_diff(i, :));
    [Q_min_temp, Q_idx_temp] = min(Q_diff(i, :));
    [I_min_temp, I_idx_temp] = min(I_diff(i, :));

    % Update arrays
    G_min_val(i) = G_min_temp; G_min_idx(i) = G_idx_temp;
    Q_min_val(i) = Q_min_temp; Q_min_idx(i) = Q_idx_temp;
    I_min_val(i) = I_min_temp; I_min_idx(i) = I_idx_temp;

  end

  %----------------------------------%
  %     Find Perturbation Vector     %
  %----------------------------------%
  % Empty arrays for perturbation size and angle
  A_perturb = zeros(length(Q_PO), 1);
  theta_perturb = zeros(length(Q_PO), 1);

  % Cycle through all Peroidic orbit values and calculate norm of distance
  % vector
  for i = 1 : length(Q_PO)
    % Get periodic orbit point
    vec_PO = [xbp_PO(i, 1), xbp_PO(i, 3)];

    % Get stable manifold point
    vec_W  = [Wq_s(Q_min_idx(i), 1), Wq_s(Q_min_idx(i), 3)];

    % Calculate difference
    vec_diff = vec_PO - vec_W;
    vec_diff = -vec_diff;
    A_diff = norm(vec_diff);

    % Angle of displacement vector
    theta_diff = mod(atan2(vec_diff(2), vec_diff(1)), 2*pi);

    % Update array
    A_perturb(i)     = A_diff;
    theta_perturb(i) = theta_diff;

  end

  %----------------%
  %     Output     %
  %----------------%
  theta_old_out     = tbp_PO;
  theta_perturb_out = theta_perturb;
  A_perturb_out     = A_perturb;

end