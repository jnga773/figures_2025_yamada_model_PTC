function run_PTC_continuation(run_new, run_old, label_old, data_PR, bcs_funcs, continuation_parameters, parameter_range, options)
  % run_PTC_continuation(run_new, run_old, label_old, data_PR, bcs_funcs, options)
  %
  % Scan through SP labels from previous run (different values of A_perturb)
  % and continue in \theta_{old} and \theta_{new}. Each run will save
  % to the 'data/run10_phase_reset_PTC/' directory.
  %
  % Parameters
  % ----------
  % run_new : string
  %     The new run identifier for the main continuation problem.
  % run_old : string
  %     The old run identifier for the sub continuation problem.
  % label_old : integer
  %     The label identifier for the previous continuation problem.
  % data_PR : struct
  %     Data structure containing the initial conditions for the trajectory segments.
  % bcs_funcs : list of functions
  %     Structure containing boundary condition functions.
  % continuation_parameters : cell
  %     Cell array containing additional parameters for the continuation.
  % parameter_range : cell
  %     Cell array containing the ranges for the continuation parameters.
  % SP_parameter : array
  %     Parameter to save SP solutions for
  % SP_values : array
  %     Array of values of theta_old to save solutions with the 'SP' label.
  % TOL : double
  %     Tolerance for the continuation. (Default: 1e-6)
  % h_min : double
  %     Minimum step size for the continuation. (default: 1e-2)
  % h0 : double
  %     Initial step size for the continuation. (default: 1e-1)
  % h_max : double
  %     Maximum step size for the continuation. (default: 1e0)
  % NAdapt : integer
  %     Number of adaptive mesh refinements. (default: 10)
  % PtMX : integer
  %     Maximum number of points in the continuation. (default: 750)
  % MaxRes : integer
  %     Maximum number of residuals allowed. (default: 10)
  % al_max : integer
  %     Maximum number of continuation steps. (default: 25)
  % NPR : integer
  %     Frequency of saved solutions. (default: 10)
  %
  % See Also
  % --------
  % coco_prob, coco_set, ode_coll2coll, apply_PR_boundary_conditions, coco_add_event, coco

  %-------------------%
  %     Arguments     %
  %-------------------%
  arguments
    run_new
    run_old
    label_old double
    data_PR struct
    bcs_funcs struct
    continuation_parameters cell = {'theta_old', 'theta_new', 'eta', 'mu_s', 'T'};
    parameter_range cell = {[0.0, 2.0], [], [-1e-4, 1e-2], [0.99, 1.01], []};

    % Optional arguments
    options.SP_parameter string = ''
    options.SP_values double = []

    % COCO Settings
    options.TOL    = 1e-6
    options.h_min  = 1e-2;
    options.h0     = 1e-1;
    options.h_max  = 1e0;
    options.NAdapt = 10;
    options.PtMX   = 750;
    options.MaxRes = 10;
    options.al_max = 25;
    options.NPR    = 10;
  end

  %----------------------------%
  %     Setup Continuation     %
  %----------------------------%
  % Set up the COCO problem
  prob = coco_prob();

  % Set tolerance
  prob = coco_set(prob, 'corr', 'TOL', options.TOL);

  % Set step sizes
  prob = coco_set(prob, 'cont', 'h_min', options.h_min);
  prob = coco_set(prob, 'cont', 'h0', options.h0);
  prob = coco_set(prob, 'cont', 'h_max', options.h_max);

  % Set adaptive meshR
  prob = coco_set(prob, 'cont', 'NAdapt', options.NAdapt);

  % Set number of steps
  prob = coco_set(prob, 'cont', 'PtMX', options.PtMX);

  % Set norm to int
  prob = coco_set(prob, 'cont', 'norm', inf);

  % Set MaxRes and al_max
  prob = coco_set(prob, 'cont', 'MaxRes', options.MaxRes);
  prob = coco_set(prob, 'cont', 'al_max', options.al_max);

  % Set frequency of saved solutions
  prob = coco_set(prob, 'cont', 'NPR', options.NPR);

  %-------------------------------------------%
  %     Continue from Trajectory Segments     %
  %-------------------------------------------%
  % Segment 1
  prob = ode_coll2coll(prob, 'seg1', run_old, label_old);
  % Segment 2
  prob = ode_coll2coll(prob, 'seg2', run_old, label_old);
  % Segment 3
  prob = ode_coll2coll(prob, 'seg3', run_old, label_old);
  % Segment 4
  prob = ode_coll2coll(prob, 'seg4', run_old, label_old);

  %------------------------------------------------%
  %     Apply Boundary Conditions and Settings     %
  %------------------------------------------------%
  % Apply all boundary conditions, glue parameters together, and
  % all that other good COCO stuff. Looking the function file
  % if you need to know more ;)
  prob = apply_boundary_conditions_PR(prob, data_PR, bcs_funcs);

  %-------------------------%
  %     Add COCO Events     %
  %-------------------------%
  if ~isempty(options.SP_values)
    prob = coco_add_event(prob, 'SP', options.SP_parameter, options.SP_values);
  end

  %-------------------------%
  %     Add COCO Events     %
  %-------------------------%
  % Run COCO continuation
  % prange = {[0.0, 2.0], [], [-1e-4, 1e-2], [0.99, 1.01], []};
  % coco(prob, run_new, [], 1, {'theta_old', 'theta_new', 'eta', 'mu_s', 'T', 'A_perturb'}, prange);
  coco(prob, run_new, [], 1, continuation_parameters, parameter_range);

end
