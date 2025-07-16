function save_fig2_data(run_in, label_in, filename_in)
  % save_fig2_data(run_in, label_in, filename_in)
  %
  % Reads periodic orbit solution data from COCO solution, calculates the
  % one-dimensional stable manifold of the "central" saddle point 'q', and
  % saves the data to filename_in.

  %-----------------------------------%
  %     Read Data: Periodic Orbit     %
  %-----------------------------------%
  % Read COCO solution
  [sol_PO, data_PO] = coll_read_solution('PO_stable.po.orb', run_in, label_in);
  
  % State space solution
  xbp_PO = sol_PO.xbp;
  % Temporal solution
  tbp_PO = sol_PO.tbp;
  % Period
  T_PO   = sol_PO.T;

  % Parameters
  p      = sol_PO.p;
  pnames = data_PO.pnames;

  %--------------------------------------%
  %     Read Data: Stationary Points     %
  %--------------------------------------%
  % Read COCO solutions
  [sol_0, ~]   = ep_read_solution('x0', run_in, label_in);
  [sol_pos, ~] = ep_read_solution('xpos', run_in, label_in);
  [sol_neg, ~] = ep_read_solution('xneg', run_in, label_in);

  % State space solutions
  x0   = sol_0.x;
  xpos = sol_pos.x;
  xneg = sol_neg.x;

  %-----------------------------%
  %     Read Data: Manifold     %
  %-----------------------------%
  % Read stable manifold solutions
  [sol1, ~] = coll_read_solution('W1', run_in, label_in);
  [sol2, ~] = coll_read_solution('W2', run_in, label_in);

  % State space solutions
  W1 = sol1.xbp;
  W2 = sol2.xbp;

  % Append to single array
  W_out = [W1; flip(W2)];

  %----------------%
  %     Output     %
  %----------------%
  % Periodic orbit solution
  data_out.xbp_PO = xbp_PO;
  data_out.tbp_PO = tbp_PO;
  data_out.T_PO   = T_PO;

  % Stationary solutions
  data_out.x0     = x0;
  data_out.xpos   = xpos;
  data_out.xneg   = xneg;

  % Stable manifold of q
  data_out.Wq_s   = W_out;

  % Parameters
  data_out.p      = p;
  data_out.pnames = pnames;

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_in, '-struct', 'data_out');

end