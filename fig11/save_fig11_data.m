function save_fig11_data(run_names_in, filename_in)
  % save_fig11_data(run_names_in, filename_in)
  %
  % Reads intersection data and saves to filename_in.

  %-------------------%
  %     Run Names     %
  %-------------------%
  runq = run_names_in.PR_move_theta;
  run1 = run_names_in.follow_theta_Wsq_part1;
  run2 = run_names_in.follow_theta_Wsq_part2;

  %--------------------------------------%
  %     Read Data: Equilibrium Point     %
  %--------------------------------------%
  % Read bifurcation data
  bd_readq = coco_bd_read(runq);

  % Get labels
  qlabs = coco_bd_labs(bd_readq, 'q_PT');

  % Read theta_old
  theta_oldq = coco_bd_val(bd_readq, 'q_PT', 'theta')';

  % Read xpos solution
  [solq, ~] = ep_read_solution('xpos', runq, qlabs(1));
  xpos = solq.x';

  % Read seg3 solution
  [sol2_1, ~] = coll_read_solution('seg2', runq, qlabs(1));
  [sol2_2, ~] = coll_read_solution('seg2', runq, qlabs(2));

  % Read state space solutions
  xbp1 = sol2_1.xbp(1, :);
  xbp2 = sol2_2.xbp(1, :);

  % Calculate distances
  vec1 = xpos - xbp1;
  vec2 = xpos - xbp2;

  % Displacement amplitude
  A1 = norm(vec1);
  A2 = norm(vec2);

  % Angle of displacement vector
  T1 = mod(atan2(vec1(3), vec1(1)), 2*pi);
  T2 = mod(atan2(vec2(3), vec2(1)), 2*pi);
  
  % Normalise
  T1 = T1 / (2 * pi);
  T2 = T2 / (2 * pi);

  % Make into arrays
  A_perturbq    = [A1; A2];
  theta_perturbq = [T1; T2];

  %-------------------------------------------%
  %     Read Data: Manifold Intersections     %
  %-------------------------------------------%
  % Read bifurcation data for both runs
  bd_read1 = coco_bd_read(run1);
  bd_read2 = coco_bd_read(run2);

  % Read theta_old values
  theta_old1     = coco_bd_col(bd_read1, 'theta')';
  theta_old2     = coco_bd_col(bd_read2, 'theta')';
  
  % Read A_perturb data
  A_perturb1     = coco_bd_col(bd_read1, 'A_perturb')';
  A_perturb2     = coco_bd_col(bd_read2, 'A_perturb')';
  
  % Read theta_perturb data
  theta_perturb1 = coco_bd_col(bd_read1, 'theta_perturb')';
  theta_perturb2 = coco_bd_col(bd_read2, 'theta_perturb')';

  % Read \Gamma coordinate date
  G_gamma1 = coco_bd_col(bd_read1, 'G_theta')';
  G_gamma2 = coco_bd_col(bd_read2, 'G_theta')';
  Q_gamma1 = coco_bd_col(bd_read1, 'Q_theta')';
  Q_gamma2 = coco_bd_col(bd_read2, 'Q_theta')';
  I_gamma1 = coco_bd_col(bd_read1, 'I_theta')';
  I_gamma2 = coco_bd_col(bd_read2, 'I_theta')';

  % Read W^{s}(q) coordinate date
  G_Wsq1 = coco_bd_col(bd_read1, 'G_Wsq')';
  G_Wsq2 = coco_bd_col(bd_read2, 'G_Wsq')';
  Q_Wsq1 = coco_bd_col(bd_read1, 'Q_Wsq')';
  Q_Wsq2 = coco_bd_col(bd_read2, 'Q_Wsq')';
  I_Wsq1 = coco_bd_col(bd_read1, 'I_Wsq')';
  I_Wsq2 = coco_bd_col(bd_read2, 'I_Wsq')';

  %----------------------%
  %     Combine Data     %
  %----------------------%
  % Combine data arrays
  % theta_old     = [theta_old1; theta_old2; theta_oldq];
  % A_perturb     = [A_perturb1; A_perturb2; A_perturbq];
  % theta_perturb = [theta_perturb1; theta_perturb2; theta_perturbq];
  theta_old     = [theta_old1; theta_old2];
  A_perturb     = [A_perturb1; A_perturb2];
  theta_perturb = [theta_perturb1; theta_perturb2];
  G_gamma       = [G_gamma1; G_gamma2];
  Q_gamma       = [Q_gamma1; Q_gamma2];
  I_gamma       = [I_gamma1; I_gamma2];
  G_Wsq         = [G_Wsq1; G_Wsq2];
  Q_Wsq         = [Q_Wsq1; Q_Wsq2];
  I_Wsq         = [I_Wsq1; I_Wsq2];

  % Mod the data
  theta_perturb = mod(theta_perturb, 1.0);

  %-----------------------------------------%
  %     Sort Data: Increasing theta_old     %
  %-----------------------------------------%
  % % Sort
  % [~, sort_idx] = sort(theta_old);
  % 
  % % Sort by theta_old
  % theta_old     = theta_old(sort_idx);
  % theta_perturb = theta_perturb(sort_idx);
  % A_perturb     = A_perturb(sort_idx);
  
  % Sort
  [~, sort_idx] = sort(theta_perturb);

  % Sort by theta_old
  theta_old     = theta_old(sort_idx);
  theta_perturb = theta_perturb(sort_idx);
  A_perturb     = A_perturb(sort_idx);
  G_gamma       = G_gamma(sort_idx);
  Q_gamma       = Q_gamma(sort_idx);
  I_gamma       = I_gamma(sort_idx);
  G_Wsq         = G_Wsq(sort_idx);
  Q_Wsq         = Q_Wsq(sort_idx);
  I_Wsq         = I_Wsq(sort_idx);

  % Get indices for 0.0 < = theta_perturb <= 0.25
  mask = (theta_perturb >= 0.0) & (theta_perturb <= 0.25);

  G_gamma = G_gamma(mask);
  Q_gamma = Q_gamma(mask);
  I_gamma = I_gamma(mask);
  G_Wsq   = G_Wsq(mask);
  Q_Wsq   = Q_Wsq(mask);
  I_Wsq   = I_Wsq(mask);

  %---------------------------------------------%
  %     Read Data: Intersection Coordinates     %
  %---------------------------------------------%

  %-------------------%
  %     Save Data     %
  %-------------------%  
  % Manifold intersection data
  data_out.theta_old     = theta_old;
  data_out.A_perturb     = A_perturb;
  data_out.theta_perturb = theta_perturb;

  data_out.xbp_gamma     = [G_gamma, Q_gamma, I_gamma];
  data_out.xbp_Wsq       = [G_Wsq, Q_Wsq, I_Wsq];

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end