%=========================================================================%
%                     YAMADA MODEL (Phase Resetting)                      %
%=========================================================================%
% We compute the phase resetting of an attracting periodic orbit of the
% Yamada model:
%                     G' = \gamma (A - G - G I) ,
%                     Q' = \gamma (B - Q - a Q I) ,
%                     I' = (G - Q - 1) I ,
% where G is the gain, Q is the absorption, and I is the intensity of the
% laser. The system is dependent on four parameters: the pump current on
% the gain, A (or A); the relative absoprtion, B and a; and the decay
% time of the gain, \gamma.
 
% Clear plots
close('all');

% Clear workspace
clear;
clc;

% Add equation/functions to path
addpath('./functions/');
% Add field functions to path
addpath('./functions/fields/');
% Add boundary condition functions to path
addpath('./functions/bcs/');
% Add SymCOCO files to path
addpath('./functions/symcoco/');

% Add continuation scripts
addpath('./continuation_scripts/');

%--------------------%
%     Parameters     %
%--------------------%
% Because we will only be looking at the (A, \gamma) plane, we will be
% setting values for a and B.
B = 5.8;
a = 1.8;

% Parameters for the periodic orbit
gamma_PO = 3.5e-2;
A_PO     = 7.4;

%-----------------------%
%     Problem Setup     %
%-----------------------%
% Parameter names
pnames = {'gamma', 'A', 'B', 'a'};

% Initial parameter values
p0 = [gamma_PO; A_PO; B; a];

% Initial point
x0 = [10; 10; 10];

% State dimensions
pdim = length(p0);
xdim = length(x0);

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Vector field: Functions
% funcs.field = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs.field = yamada_symbolic();

% Adjoint equations: Functions (for floquet_mu and floquet_wnorm)
% funcs.VAR = {@VAR};
funcs.VAR = VAR_symbolic();

% Phase Reset Segment 1: Functions
% func.seg1 = {@func_seg1};
funcs.seg1 = func_seg1_symbolic();

% Phase Reset: Segment 2
% funcs.seg2 = {@func_seg2};
funcs.seg2 = func_seg2_symbolic();

% Phase Reset: Segment 3
% funcs.seg3 = {@func_seg3};
funcs.seg3 = func_seg3_symbolic();

% Phase Reset: Segment 4
% funcs.seg4 = {@func_seg4};
funcs.seg4 = func_seg4_symbolic();

% Boundary conditions: Period
% bcs_funcs.bcs_T = {@bcs_T};
bcs_funcs.bcs_T = bcs_T_symbolic();

% Boundary conditions: Periodic orbit
% bcs_funcs.bcs_PO = {@bcs_PO};
bcs_funcs.bcs_PO = bcs_PO_symbolic();

% Boundary conditions: Floquet multipliers
% bcs_funcs.bcs_VAR = {@bcs_VAR};
bcs_funcs.bcs_VAR = bcs_VAR_symbolic();

% Boundary conditions: Phase-resetting segments
% bcs_funcs.bcs_PR = {@bcs_PR};
bcs_funcs.bcs_PR = bcs_PR_symbolic();

%=========================================================================%
%%                   CALCULATE INITIAL PERIODIC ORBIT                    %%
%=========================================================================%
% Using ODE45, we compute a guess solution to a stable periodic orbit. We
% then feed this as an initial solution to the 'PO' toolbox. Finally, we
% "rotate" the head-point and use this to confirm a solution of a periodic
% orbit, where the first point corresponds to max(G).

%-------------------------------------------------------------------------%
%%                 Confirm ODE45 Periodic Orbit Solution                 %%
%-------------------------------------------------------------------------%
% Calculate the periodic orbit using MATLAB's ode45 function.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.ode45_PO = 'run01_initial_PO_ode45';
run_new = run_names.ode45_PO;

% Print to console
fprintf("~~~ Initial Periodic Orbit: First Run ~~~ \n");
fprintf('Find new periodic orbit \n');
fprintf('Run name: %s \n', run_new);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_ode45 = calc_initial_solution_ODE45(x0, p0, funcs.field);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set initial guess to 'coll'
prob = ode_isol2po(prob, '', funcs.field{:}, ...
                   data_ode45.t, data_ode45.x, pnames, p0);

% Add equilibrium points for non trivial steady states
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, ...
                   data_ode45.xpos, p0);
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, ...
                   data_ode45.xneg, p0);
prob = ode_isol2ep(prob, 'x0', funcs.field{:}, ...
                   data_ode45.x0, p0);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Glue parameters and apply boundary condition
prob = glue_parameters_PO(prob);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
prob = coco_add_event(prob, 'PO_PT', 'A', A_PO);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'});

%-------------------------------------------------------------------------%
%%                   Re-Solve for Rotated Perioid Orbit                  %%
%-------------------------------------------------------------------------%
% Using previous parameters and MATLAB's ode45 function, we solve for an
% initial solution to be fed in as a periodic orbit solution.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.initial_PO = 'run02_initial_periodic_orbit';
run_new = run_names.initial_PO;
% Which run this continuation continues from
run_old = run_names.ode45_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');
label_old = label_old(1);

% Print to console
fprintf("~~~ Initial Periodic Orbit: Sixth Run ~~~ \n");
fprintf('Find new periodic orbit \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Calculate Solution     %
%----------------------------%
% Calculate dem tings
data_soln = calc_initial_solution_PO(run_old, label_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set NTST mesh 
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set PtMX steps
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Set initial guess to 'coll'
prob = ode_isol2coll(prob, 'initial_PO', funcs.field{:}, ...
                     data_soln.t, data_soln.x, pnames, data_soln.p);

% Add equilibrium points for non trivial steady states
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
prob = ode_ep2ep(prob, 'x0',   run_old, label_old);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Glue parameters and apply boundary condition
prob = apply_boundary_conditions_PO(prob, bcs_funcs.bcs_PO);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Event for A = 7.5
prob = coco_add_event(prob, 'PO_PT', 'A', data_soln.p(2));

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'A', 'gamma'});

%=========================================================================%
%%            Compute Floquet Bundle at Zero Phase Point (mu)            %%
%=========================================================================%
% We now add the adjoint function and Floquet boundary conditions to
% compute the adjoint (left or right idk) eigenvectors and eigenvalues.
% This will give us the perpendicular vector to the tangent of the periodic
% orbit. However, this will only be for the eigenvector corresponding to
% the eigenvalue \mu = 1. Hence, here we continue in \mu (mu_s) until
% mu_s = 1.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.compute_floquet_1 = 'run03_compute_floquet_bundle_1_mu';
run_new = run_names.compute_floquet_1;
% Which run this continuation continues from
run_old = run_names.initial_PO;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');

% Print to console
fprintf("~~~ Floquet Bundle: First Run ~~~ \n");
fprintf('Calculate Floquet bundle (mu) \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%--------------------------%
%     Calculate Things     %
%--------------------------%
data_adjoint = calc_initial_solution_VAR(run_old, label_old);

%------------------------------------%
%     Setup Floquet Continuation     %
%------------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2);
prob = coco_set(prob, 'cont', 'h0', 1e-2);
prob = coco_set(prob, 'cont', 'h_max', 1e-2);

% Set PtMX
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set NTST
prob = coco_set(prob, 'coll', 'NTST', 50);

% Set NAdapt
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Add segment as initial solution
prob = ode_isol2coll(prob, 'adjoint', funcs.VAR{:}, ...
                     data_adjoint.t0, data_adjoint.x0, ...
                     data_adjoint.pnames, data_adjoint.p0);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply boundary conditions
prob = apply_boundary_conditions_VAR(prob, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Add event
prob = coco_add_event(prob, 'mu=1', 'mu_s', 1.0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'mu_s', 'w_norm', 'T'} , [0.0, 1.1]);

%-------------------------------------------------------------------------%
%%          Compute Floquet Bundle at Zero Phase Point (w_norm)          %%
%-------------------------------------------------------------------------%
% Having found the solution (branching point 'BP') corresponding to
% \mu = 1, we can continue in the norm of the vector w (w_norm), until the
% norm is equal to zero. Then we will have the correct perpendicular
% vector.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.compute_floquet_2 = 'run04_compute_floquet_bundle_2_w';
run_new = run_names.compute_floquet_2;
% Which run this continuation continues from
run_old = run_names.compute_floquet_1;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'BP');
label_old = label_old(1);

% Print to console
fprintf("~~~ Floquet Bundle: Second Run ~~~ \n");
fprintf('Calculate Floquet bundle (w_norm) \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%------------------------------------%
%     Setup Floquet Continuation     %
%------------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set number of PtMX steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set number of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Continue coll from previous branching point
% prob = ode_BP2coll(prob, 'adjoint', run_old, label_old);
prob = coco_set(prob, 'cont', 'branch', 'switch');
prob = ode_coll2coll(prob, 'adjoint', run_old, label_old);

%------------------------------------------------%
%     Apply Boundary Conditions and Settings     %
%------------------------------------------------%
% Apply boundary conditions
prob = apply_boundary_conditions_VAR(prob, bcs_funcs);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Add event when w_norm = 1
prob = coco_add_event(prob, 'NORM1', 'w_norm', 1.0);

%------------------%
%     Run COCO     %
%------------------%
% Run COCO continuation
coco(prob, run_new, [], 1, {'w_norm', 'mu_s', 'T'}, {[-1e-4, 1.1], [], []});

%=========================================================================%
%%                   CALCULATE PHASE RESET SOLUTIONS                     %%
%=========================================================================%
% We compute a set of phase reset problems. We first continue in the
% parameter 'A_perturb', that is, the amplitude of the perturbation.
% After this, we then compute some phase transition curves (PTCs), by
% continuing in 'theta_old' and 'theta_new'.

%-------------------------------------------------------------------------%
%%                   Increasing Pertubation Amplitude                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.phase_reset_orbit = 'run05_phase_reset_orbit';
run_new = run_names.phase_reset_orbit;
% Which run this continuation continues from
run_old = run_names.compute_floquet_2;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'NORM1');
label_old = label_old(1);

% Print to console
fprintf("~~~ Phase Transition Curve: First Run ~~~ \n");
fprintf('Move along periodic orbit \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%-------------------%
%     Read Data     %
%-------------------%
% Set periodicity
k = 30;

% Set perturbation direction to be d = (1, 0, 1) / sqrt(2)
theta_perturb = deg2rad(135);

% Set initial conditions from previous solutions
data_PR = calc_initial_solution_PR(run_old, label_old, k, theta_perturb);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 5e-5);
prob = coco_set(prob, 'cont', 'h0', 1e-3);
prob = coco_set(prob, 'cont', 'h_max', 1e0);

% Set adaptive mesh
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set number of steps
prob = coco_set(prob, 'cont', 'PtMX', 100);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Set norm to int
prob = coco_set(prob, 'cont', 'norm', inf);

% Set MaxRes and al_max
prob = coco_set(prob, 'cont', 'MaxRes', 10);
prob = coco_set(prob, 'cont', 'al_max', 25);

%------------------%
%     Set NTST     %
%------------------%
% In calc_PR_initial conditions, we define segment 4 as having 'k' periods,
% where 'k' is an integer. This is the perturbed segment, that may have to
% orbit the unperturbed periodic orbit many times before "resetting". Hence
% we have set the NTST for this segment (NTST(4)) as k * 50.
NTST(1) = 30;
NTST(2) = 10;
NTST(3) = 10;
NTST(4) = 30 * k;

prob = coco_set(prob, 'seg1.coll', 'NTST', NTST(1));
prob = coco_set(prob, 'seg2.coll', 'NTST', NTST(2));
prob = coco_set(prob, 'seg3.coll', 'NTST', NTST(3));
prob = coco_set(prob, 'seg4.coll', 'NTST', NTST(4));

%------------------------------------%
%     Create Trajectory Segments     %
%------------------------------------%
% Here we set up the four segments to calculate the phase resetting curve:
% Segment 1 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point where the perturbed trajectory 
%             comes close to the periodic orbit (at theta_new).
prob = ode_isol2coll(prob, 'seg1', funcs.seg1{:}, ...
                     data_PR.t_seg1, data_PR.x_seg1, data_PR.p0);

% Segment 2 - Trajectory segment of the periodic from the end of Segment 1
%             (at theta_new) back to the zero-phase point (gamma_0).
prob = ode_isol2coll(prob, 'seg2', funcs.seg2{:}, ...
                     data_PR.t_seg2, data_PR.x_seg2, data_PR.p0);

% Segment 3 - Trajectory segment of the periodic orbit from the zero-phase
%             point (gamma_0) to the point at which the perturbation is
%             applied (theta_old).
prob = ode_isol2coll(prob, 'seg3', funcs.seg3{:}, ...
                     data_PR.t_seg3, data_PR.x_seg3, data_PR.p0);   

% Segment 4 - Trajectory segment of the perturbed trajectory, from
%             theta_old to theta_new.
prob = ode_isol2coll(prob, 'seg4', funcs.seg4{:}, ...
                     data_PR.t_seg4, data_PR.x_seg4, data_PR.p0);

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
% Save solution at phase along \Gamma where there WILL BE an intersection
% with the stable manifold of q.
prob = coco_add_event(prob, 'SP', 'theta_old', [0.339279, 1.339279]);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and parameter range
pcont = {'theta_old', 'theta_new', ...
         'eta', 'mu_s', 'T'};
prange = {[1.0, 2.0], [1.0, 2.0], ...
          [-1e-4, 1e-2], [0.99, 1.01], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                   Increasing Pertubation Amplitude                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.phase_reset_perturbation = 'run06_phase_reset_perturbation';
run_new = run_names.phase_reset_perturbation;
% Which run this continuation continues from
run_old = run_names.phase_reset_orbit;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SP');
label_old = label_old(1);

% Print to console
fprintf("~~~ Phase Transition Curve: Second Run ~~~ \n");
fprintf('Increase perturbation amplitude \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from point %d in run: %s \n', label_old, run_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set tolerance
% prob = coco_set(prob, 'corr', 'TOL', 5e-7);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-1);
prob = coco_set(prob, 'cont', 'h0', 1e0);
prob = coco_set(prob, 'cont', 'h_max', 1e1);

% Set adaptive meshR
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set number of steps
prob = coco_set(prob, 'cont', 'PtMX', 1000);

% Set norm to int
prob = coco_set(prob, 'cont', 'norm', inf);

% % Set MaxRes and al_max
% prob = coco_set(prob, 'cont', 'MaxRes', 10);
% prob = coco_set(prob, 'cont', 'al_max', 25);

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
% Read intensity for intersection with I=0 plane
I_intsct = coco_bd_val(coco_bd_read(run_old), label_old, 'I_theta_n');

% List of perturbation amplitudes to save solutions for
SP_values = [0.1, 0.724587, 1.5257, 4.0];
prob = coco_add_event(prob, 'SP', 'A_perturb', SP_values);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Set continuation parameters and parameter range
pcont = {'A_perturb', 'theta_new', ...
         'eta', 'mu_s', 'T'};
prange = {[0.0, 4.0], [], ...
          [-1e-4, 1e-2], [0.99, 1.01], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%             Directional Transition Curve (DTC) - Multiple             %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.phase_transition_curve = 'run07_phase_reset_DTC_scan';
run_new = run_names.phase_transition_curve;
% Which run this continuation continues from
run_old = run_names.phase_reset_perturbation;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SP');

% Print to console
fprintf("~~~ Phase Reset: Second Run ~~~ \n");
fprintf('Calculate phase transition curve \n');
fprintf('Run name: %s \n', run_new);
fprintf('Continuing from SP points in run: %s \n', run_old);

%---------------------------------%
%     Cycle through SP labels     %
%---------------------------------%
% Set number of threads
M = 0;
% parfor (run = 1 : length(label_old), M)
for run = 1
  % Label for this run
  this_run_label = label_old(run);

  % Data directory for this run
  fprintf('\n Continuing from point %d in run: %s \n', this_run_label, run_old);

  this_run_name = {run_new; sprintf('run_%02d', run)};

  % Saved solution points for theta_old
  SP_values = [0.5, 1.5];

  % Continuation parameters
  continuation_parameters = {'theta_perturb', 'theta_new', 'eta', 'mu_s', 'T'};
  % Parameter range for continuation
  parameter_range = {[0.0, 2*pi], [], [-1e-4, 1e-2], [0.99, 1.01], []};

  % Run continuation
  run_PTC_continuation(this_run_name, run_old, this_run_label, data_PR, bcs_funcs, ...
                       continuation_parameters, parameter_range, ...
                       SP_parameter='theta_new', SP_values=SP_values, ...
                       h_min=1e-3, h0=5e-1, h_max=1e1, ...
                       PtMX=2000, NPR=50, NAdapt=20);

end

%=========================================================================%
%%                          SAVE AND PLOT DATA                           %%
%=========================================================================%
%-------------------%
%     Save Data     %
%-------------------%
% Save data for Figure 12
save_fig12_data(run_new, '../data_files/fig12_data.mat');

%----------------------%
%     Plot Figures     %
%----------------------%
% Run plotting scripts
plot_fig12;

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%