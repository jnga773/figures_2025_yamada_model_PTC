function data_out = insert_large_time_segment(run_in, label_in, funcs_in)
  % min_idx = insert_large_time_segment(run_in, label_in, funcs_in)
  %
  % Reads the periodic orbit solution from solution [label_old] of
  % [run_old], and finds the segment of the state-space solution
  % closest to the equilibrium point.
  %
  % With this point found, we insert a large time segment to
  % "trick" the periodic orbit into having a larger period.
  %
  % Parameters
  % ----------
  % run_in: string
  %     The string identifier for the previous COCO run that we will
  %     read information from.
  % label_in: int
  %     The solution label we will read the data from
  % funcs_in: cell
  %     Cell of field functions and Jacobians etc.
  %
  % Returns
  % -------
  % min_idx : data structure
  %     Contains the state space solution, temporal solution and
  %     parameters.
  %
  % See Also
  % --------
  % coll_read_solution

  % Yamada function
  func_yamada = funcs_in{1};

  %-------------------%
  %     Read Data     %
  %-------------------%
  % Read solution with maximum period
  [sol, data] = coll_read_solution('po.orb', run_in, label_in);

  % State-space solution
  xbp_PO = sol.xbp;
  tbp_PO = sol.tbp;

  % Period
  T_PO   = sol.T;

  % Parameters
  p0     = sol.p;
  % Parameter names
  pnames = data.pnames;

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Evaluate vector field at basepoints
  f = func_yamada(xbp_PO', repmat(p0, [1, size(xbp_PO, 1)]));

  % Extract the discretisation points corresponding to the minimum value of
  % the norm of the vector field along the longest-period periodic orbit.
  % Find basepoint closest to equilibrium
  f_norm = sqrt(sum(f .* f, 1));
  f_norm = f_norm';

  % Find min and index
  [~, min_idx] = min(f_norm);

  % Then insert a time segment that is a large multiple of the orbit
  % period immediately following the discretisation point.
  scale = 1000;

  % fprintf('Maximum Period from run ''%s'', T = %f \n', run_new, T);
  % fprintf('Scaled period is T'' = %d x %f = %f \n', scale, T, scale * T);

  % Crank up period by factor scale
  t_soln = [tbp_PO(1:min_idx,1);
            T_PO * (scale - 1) + tbp_PO(min_idx+1:end,1)];

  % Approximate equilibrium point
  x0_approx = sol.xbp(min_idx, :);

  %----------------%
  %     Output     %
  %----------------%
  % Output data structure
  data_out.p          = p0;
  data_out.pnames     = pnames;
  data_out.NTST       = data.coll.NTST;

  data_out.tbp_extend = t_soln;
  data_out.xbp_extend = xbp_PO;

  data_out.x0         = x0_approx;

end