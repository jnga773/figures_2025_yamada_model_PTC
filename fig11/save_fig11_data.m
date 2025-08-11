function save_fig11_data(run_names_in, filename_in)
  % save_fig11_data(run_names_in, filename_in)
  %
  % Reads intersection data and saves to filename_in.
  %
  % Parameters
  % ----------
  % run_names_in : struct
  %    List of all run name variables.
  % filename_in : char
  %    Name of output data file to save data to.
  %
  % See Also
  % --------
  % coco_bd_read, coco_bd_labs, coll_read_solution, ep_read_solution

  %-------------------%
  %     Run Names     %
  %-------------------%
  % Periodic orbit run name
  run_PO = run_names_in.initial_PO_COLL;
  % Intersection point run names
  runq   = run_names_in.PR_move_theta;
  run1   = run_names_in.follow_theta_Wsq_part1;
  run2   = run_names_in.follow_theta_Wsq_part2;

  %------------------------------------------%
  %     Read Initial Periodic Orbit Data     %
  %------------------------------------------%
  % Bifurcation data
  bd_PO = coco_bd_read(run_PO);

  % Get solution label
  label_PO = coco_bd_labs(bd_PO, 'PO_PT');

  % Read 'initial_PO' COLL data
  [sol_PO, data_PO] = coll_read_solution('initial_PO', run_PO, label_PO);

  % State space solution
  xbp_PO = sol_PO.xbp;
  % Temporal solution
  tbp_PO = sol_PO.tbp;
  % Period
  T_PO   = sol_PO.T;

  % Parameters
  p      = sol_PO.p;
  pnames = data_PO.pnames;

  %-------------------------------------%
  %     Read Data: Stationary Point     %
  %-------------------------------------%
  % Read 'xpos' EP data
  [sol_pos, ~] = ep_read_solution('xpos', run_PO, label_PO);

  % Stationary point
  xpos = sol_pos.x;

  %--------------------------------------------%
  %     Calculate Stable Manifold for Plot     %
  %--------------------------------------------%
  % Calculate that thang
  Wq_s = calc_stable_manifold(run_PO, label_PO);

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
  labs1 = coco_bd_labs(bd_read1, 'SP');
  labs2 = coco_bd_labs(bd_read2, 'SP');

  % Intersection points
  theta_old_SP     = [coco_bd_val(bd_read1, labs1, 'theta'), coco_bd_val(bd_read2, labs2, 'theta')]';
  A_perturb_SP     = [coco_bd_val(bd_read1, labs1, 'A_perturb'), coco_bd_val(bd_read2, labs2, 'A_perturb')]';
  theta_perturb_SP = [coco_bd_val(bd_read1, labs1, 'theta_perturb'), coco_bd_val(bd_read2, labs2, 'theta_perturb')]';
  I_theta_SP       = [coco_bd_val(bd_read1, labs1, 'I_theta'), coco_bd_val(bd_read2, labs2, 'I_theta')]';

  % Get \Gamma and W^{s}(q) coordinates
  G_theta_SP = [coco_bd_val(bd_read1, labs1, 'G_theta'), coco_bd_val(bd_read2, labs2, 'G_theta')]';
  Q_theta_SP = [coco_bd_val(bd_read1, labs1, 'Q_theta'), coco_bd_val(bd_read2, labs2, 'Q_theta')]';
  I_theta_SP = [coco_bd_val(bd_read1, labs1, 'I_theta'), coco_bd_val(bd_read2, labs2, 'I_theta')]';
  G_Wsq_SP   = [coco_bd_val(bd_read1, labs1, 'G_Wsq'), coco_bd_val(bd_read2, labs2, 'G_Wsq')]';
  Q_Wsq_SP   = [coco_bd_val(bd_read1, labs1, 'Q_Wsq'), coco_bd_val(bd_read2, labs2, 'Q_Wsq')]';
  I_Wsq_SP   = [coco_bd_val(bd_read1, labs1, 'I_Wsq'), coco_bd_val(bd_read2, labs2, 'I_Wsq')]';

  % Turn into data arrays
  xbp_gamma_SP = [G_theta_SP, Q_theta_SP, I_theta_SP];
  xbp_Wsq_SP   = [G_Wsq_SP, Q_Wsq_SP, I_Wsq_SP];

  %-------------------%
  %     Save Data     %
  %-------------------%
  % Parameters
  data_out.p                = p;
  data_out.pnames           = pnames;

  % Periodic orbit data
  data_out.xbp_PO           = xbp_PO;
  data_out.tbp_PO           = tbp_PO;
  data_out.T_PO             = T_PO;

  % Stationary point
  data_out.xpos             = xpos;
  data_out.Wq_s             = Wq_s;

  % Manifold intersection data
  data_out.theta_old        = theta_old;
  data_out.A_perturb        = A_perturb;
  data_out.theta_perturb    = theta_perturb;
  
  % Highlighted points along \Gamma and W^{s}(q)
  data_out.xbp_gamma        = [G_gamma, Q_gamma, I_gamma];
  data_out.xbp_Wsq          = [G_Wsq, Q_Wsq, I_Wsq];

  % Special intersection points
  data_out.theta_old_SP     = theta_old_SP;
  data_out.A_perturb_SP     = A_perturb_SP;
  data_out.theta_perturb_SP = theta_perturb_SP;
  data_out.I_theta_SP       = I_theta_SP;

  data_out.xbp_gamma_SP     = xbp_gamma_SP;
  data_out.xbp_Wsq_SP       = xbp_Wsq_SP;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end

function x_out = calc_stable_manifold(run_in, label_in)
  % x_out = calc_stable_manifold(run_in, label_in)
  %
  % Calculate the stable manifold of the equilibrium point in the middle of
  % the periodic orbit.
  %
  % Parameters
  % ----------
  % run_in : char
  %     The run identifier for the continuation problem.
  % label_in : float
  %     The solution label for the continuation problem.
  %
  % Returns
  % -------
  % x_out : array
  %    State-space solution of the stable manifold of q.
  %
  % See Also
  % --------
  % ep_read_solution, ode45

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_in char
    label_in double
  end

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read EP solution
  [sol, data]   = ep_read_solution('xpos', run_in, label_in);

  % Equilibrium point
  x_pos = sol.x;
  % Parameters
  p     = sol.p;

  % Function handles
  fhan    = data.fhan;
  DFDXhan = data.dfdxhan;
  
  %------------------------------%
  %     Calculate EigenStuff     %
  %------------------------------%
  % Jacobian
  J_stable = DFDXhan(x_pos, p);

  % Calculate eigenvalues and eigenvectors
  [eigvec, ~] = eig(J_stable);

  % Indices for stable eigenvectors (eigval < 0)
  % stable_index = find(diag(eigval) < 0);
  stable_index = 3;

  % Stable eigenvector
  vec_s = eigvec(:, stable_index);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps1 = -0.01;
  % Time span
  t_span1 = -16.5 : 0.01 : 0.0;
  t_span1 = flip(t_span1);

  % Initial vector
  x_init1 = x_pos + (eps1 * vec_s);

  % Integrate using ode45
  [~, W1] = ode45(@(t_in, x0_in) fhan(x0_in, p), t_span1, x_init1);

  %------------------------------------------------%
  %     Calculate Things: Positive y direction     %
  %------------------------------------------------%
  % Small distance
  eps2 = 0.01;
  % Time span
  t_span2 = -25.0 : 0.01 : 0.0;
  t_span2 = flip(t_span2);

  % Initial vector
  x_init2 = x_pos + (eps2 * vec_s);

  % Integrate using ode45
  [~, W2] = ode45(@(t_in, x0_in) fhan(x0_in, p), t_span2, x_init2);

  %----------------%
  %     Output     %
  %----------------%
  x_out = [flip(W2); W1];

end