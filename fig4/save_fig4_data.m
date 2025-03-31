function save_fig8_data(run_in, filename_in)
  % save_fig8_data(run_in, filename_in)
  %
  % Reads PTC scan data from run_in and saves data to filename_in

  %------------------------------------------%
  %     Read Initial Periodic Orbit Data     %
  %------------------------------------------%
  % Run string identifier
  run_PO = 'run02_initial_periodic_orbit';
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
  dir_data = sprintf('./data/%s/', run_in);
  % List all directories
  dirs = dir(dir_data);
  % Remove ./ and ../
  dirs = dirs(~ismember({dirs.name}, {'.', '..', '.DS_Store'}));
  % Sub folder names
  dir_sub = {dirs.name};

  % Empty cells
  A_perturb = cell(1, 2);
  sol3      = cell(1, 2);
  sol4      = cell(1, 2);

  % Cycle through data sub directories
  for i = 1 : length(dir_sub)
    % Run name
    sub_run_name = {run_in, dir_sub{i}};

    % Bifurcation data
    bd_PR = coco_bd_read(sub_run_name);

    % Get labels
    label_PR = coco_bd_labs(bd_PR, 'SP');
    label_PR = label_PR(1);

    % Get theta_old values
    theta_old = coco_bd_val(bd_PR, label_PR, 'theta_old');

    % Get A_perturb value
    A_perturb_read = coco_bd_val(bd_PR, label_PR, 'A_perturb');
  
    % Get periodicity
    k = coco_bd_val(bd_PR, label_PR, 'k');

    % Read segment 3 solution
    [sol3_read, ~] = coll_read_solution('seg3', sub_run_name, label_PR);

    % Read segment 4 solution
    [sol4_read, ~] = coll_read_solution('seg4', sub_run_name, label_PR);

    % Append to arrays
    A_perturb{i} = A_perturb_read;
    sol3{i}      = sol3_read;
    sol4{i}      = sol4_read;
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
  data_out.A_perturb_run1 = A_perturb{1};
  data_out.A_perturb_run2 = A_perturb{2};
  data_out.k              = k;
  data_out.T_PO           = T_PO;

  % Plotting data
  data_out.xbp_PO_plot    = xbp_PO_plot;
  data_out.tbp_PO_plot    = tbp_PO_plot;

  % Save to Matlab .mat data structure
  save(filename_in, '-struct', 'data_out');

end