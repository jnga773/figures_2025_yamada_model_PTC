function save_fig3_data(run_names_in, filename_in)
  % save_fig3_data(run_names_in, filename_in)
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
  run_PO = run_names_in.initial_PO_COLL;
  run_PR = run_names_in.PR_perturbation;

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

  %--------------------------------%
  %     Read Data: Phase Resets    %
  %--------------------------------%
  % Bifurcation data
  bd_PR = coco_bd_read(run_PR);

  % Get labels
  label_PR = coco_bd_labs(bd_PR, 'SP');

  % Read perturbation amplitude
  A_perturb = coco_bd_val(bd_PR, 'SP', 'A_perturb');
  % Read theta_old phase
  theta_old = coco_bd_val(bd_PR, 'SP', 'theta_old');
  theta_old = theta_old(1);
  % Read periodicity
  k         = coco_bd_val(bd_PR, 'SP', 'k');
  k         = k(1);

  % Empty arrays for data
  sol3 = cell(length(label_PR), 1);
  sol4 = cell(length(label_PR), 1);

  % Cycle through data sub directories
  for idx = 1 : length(label_PR)
    % Read segment 3 solution
    [sol3_read, ~] = coll_read_solution('seg3', run_PR, label_PR(idx));

    % Read segment 4 solution
    [sol4_read, ~] = coll_read_solution('seg4', run_PR, label_PR(idx));

    % Append to arrays
    sol3{idx} = sol3_read;
    sol4{idx} = sol4_read;
  end

  % Get segment4 state space solution
  xbp3_run1 = sol3{1}.xbp;
  xbp3_run2 = sol3{2}.xbp;
  xbp4_run1 = sol4{1}.xbp;
  xbp4_run2 = sol4{2}.xbp;

  % Get segment 4 temporal solution
  tbp4_run1 = sol4{1}.tbp;
  tbp4_run2 = sol4{2}.tbp;

  %----------------------------------------%
  %     Copy Unperturned Orbit k Times     %
  %----------------------------------------%
  % Time data
  tbp_normalised = tbp_PO / T_PO;

  % Setup plotting data
  xbp_PO_plot = xbp_PO(:, 3);
  tbp_PO_plot = tbp_normalised;

  % Copy original periodic orbit data
  for i = 2: k
    xbp_PO_plot = [xbp_PO_plot(1:end-1); xbp_PO(:, 3)];
    tbp_PO_plot = [tbp_PO_plot(1:end-1); tbp_PO_plot(end) + tbp_normalised];
  end

  % Multiply segment time solutions
  tbp4_run1 = k * tbp4_run1;
  tbp4_run2 = k * tbp4_run2;

  % Normalise by period
  tbp4_run1 = tbp4_run1 / T_PO;
  tbp4_run2 = tbp4_run2 / T_PO;

  %-------------------%
  %     Save Data     %
  %-------------------%
  % State space solutions
  data_out.xbp_PO         = xbp_PO;
  data_out.xbp4_run1      = xbp4_run1;
  data_out.xbp4_run2      = xbp4_run2;
  data_out.xbp3_run1      = xbp3_run1;
  data_out.xbp3_run2      = xbp3_run2;

  % Stationary point
  data_out.xpos           = xpos;

  % Temporal data
  data_out.tbp_PO         = tbp_PO;
  data_out.tbp4_run1      = tbp4_run1;
  data_out.tbp4_run2      = tbp4_run2;

  % Parameters
  data_out.p              = p;
  data_out.pnames         = pnames;
  data_out.theta_old      = theta_old;
  data_out.A_perturb_run1 = A_perturb(1);
  data_out.A_perturb_run2 = A_perturb(2);
  data_out.k              = k(1);
  data_out.T_PO           = T_PO;
  data_out.theta_old      = theta_old(1);

  % Plotting data
  data_out.xbp_PO_plot    = xbp_PO_plot;
  data_out.tbp_PO_plot    = tbp_PO_plot;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end