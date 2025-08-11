function save_fig6_data(run_names_in, filename_in)
  % save_fig6_data(run_names_in, filename_in)
  %
  % Reads PTC scan data from run_PO and run_PR_PTC_multi and
  % saves data to filename_in.
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
  %     Arguments     %
  %-------------------%
  arguments
    run_names_in struct
    filename_in char
  end

  %-------------------%
  %     Run Names     %
  %-------------------%
  % Periodic orbit run name
  run_PO  = run_names_in.initial_PO_COLL;
  run_PTC = run_names_in.PR_PTC_multi;

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

  %------------------------%
  %     Read Data: PTCs    %
  %------------------------%
  % Folder name
  dir_data = sprintf('./data/%s/', run_PR_PTC_multi);
  % List all directories
  dirs = dir(dir_data);
  % Remove ./ and ../
  dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
  % Sub folder names
  dir_sub = {dirs.name};

  %------------------------------------%
  %     Cycle Through Data Folders     %
  %------------------------------------%
  % Empty array for A_perturb values
  A_perturb = zeros(length(dir_sub), 1);
  theta_old = cell(1, length(dir_sub));
  theta_new = cell(1, length(dir_sub));

  for i = 1 : length(dir_sub)
  % for i = 1 : 1
    % Run name
    sub_run_name = {run_PR_PTC_multi, dir_sub{i}};

    % Bifurcation data
    bd_read = coco_bd_read(sub_run_name);

    % Read A_perturb value
    A_perturb_read = coco_bd_val(bd_read, 1, 'A_perturb');
    % Update array
    A_perturb(i) = A_perturb_read;

    % Read PTC data
    theta_old_read = coco_bd_col(bd_read, 'theta_old');
    theta_new_read = coco_bd_col(bd_read, 'theta_new');

    % Get unique values
    [~, unique_idx] = unique(theta_old_read);
    theta_old_unique = theta_old_read(unique_idx);
    theta_new_unique = theta_new_read(unique_idx);

    % Split theta_old_read into theta_old < 1 and theta_old > 1
    theta_old_lt1 = theta_old_unique(theta_old_unique <= 1.0);
    theta_old_gt1 = theta_old_unique(theta_old_unique > 1.0) - 1.0;

    % Split theta_new_read into theta_old < 1 and theta_old > 1
    theta_new_lt1 = theta_new_unique(theta_old_unique <= 1.0);
    theta_new_gt1 = theta_new_unique(theta_old_unique > 1.0);

    % Check if theta_old_lt1 and theta_old_gt1 cover entire 0 -> 1 range
    lt1_check = false;
    gt1_check = false;

    if round(min(theta_old_lt1), 3) == 0.0 && round(max(theta_old_lt1), 3) == 1.0
      lt1_check = true;
    end
    if round(min(theta_old_gt1), 3) == 0.0 && round(max(theta_old_gt1), 3) == 1.0
      gt1_check = true;
    end

    % If both are true, save the lt1 data
    if lt1_check && gt1_check
      theta_old{i} = theta_old_lt1;
      theta_new{i} = theta_new_lt1;
    end

    % If only one is true, save that data
    if lt1_check && ~gt1_check
      theta_old{i} = theta_old_lt1;
      theta_new{i} = theta_new_lt1;
    elseif ~lt1_check && gt1_check
      theta_old{i} = theta_old_gt1;
      theta_new{i} = theta_new_gt1;
    end

    % If both are false, merge the two together with a NaN in between
    if ~lt1_check && ~gt1_check
      % Append data
      theta_old{i} = [theta_old_gt1, NaN, theta_old_lt1];
      theta_new{i} = [theta_new_gt1, NaN, theta_new_lt1];
    end

    % Check if theta_new starts in the fundamental domain. If not, shift it
    if theta_old{i}(1) < 0.0
      theta_new{i} = theta_new{i} + 1.0;
    elseif theta_old{i}(1) > 1.0
      theta_new{i} = theta_new{i} - 1.0;
    end

  end

  %-------------------%
  %     Save Data     %
  %-------------------%
  % Periodic orbit data
  data_out.xbp_PO     = xbp_PO;
  data_out.tbp_PO     = tbp_PO;
  data_out.T_PO       = T_PO;

  % Stationary point
  data_out.xpos       = xpos;
  data_out.Wq_s       = Wq_s;
  
  % PTC data
  data_out.theta_old  = theta_old;
  data_out.theta_new  = theta_new;
  data_out.A_perturb  = A_perturb;

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
  % label_in : double
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