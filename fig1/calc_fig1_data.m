%=========================================================================%
%                   YAMADA MODEL (Bifurcation Diagram)                    %
%=========================================================================%
% We compute the two parameter bifurcation diagram (in A and \gamma) for
% the Yamada model:
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
addpath('../COCO_files/');
% Add field functions to path
addpath('../COCO_files/fields/');
% Add boundary condition functions to path
addpath('../COCO_files/bcs/');
% Add SymCOCO files to path
addpath('../COCO_files/symcoco/');

% Add continuation scripts
addpath('./continuation_scripts/homoclinic_approx/');
addpath('./continuation_scripts/double_limit_cycle/');
addpath('./continuation_scripts/lins_method/');

%--------------------%
%     Parameters     %
%--------------------%
% Because we will only be looking at the (A, \gamma) plane, we will be
% setting values for a and B.
B = 5.8;
a = 1.8;

% Set some initial values for \gamma and A
gamma = 0.10;
A = 6.6;

%-----------------------%
%     Problem Setup     %
%-----------------------%
% Parameter names
pnames = {'gamma', 'A', 'B', 'a'};

% Initial parameter values
p0 = [gamma; A; B; a];

% Initial state values
x0 = [A; B; 0];

% Parameter ranges
gamma_range = [0.0, 0.25];
A_range = [5.0, 11.0];

% State dimensions
pdim = length(p0);
xdim = length(x0);

%-------------------------%
%     Functions Lists     %
%-------------------------%
% Yamada model equations
% funcs.field = {@yamada, @yamada_DFDX, @yamada_DFDP};
funcs.field = yamada_symbolic();

% Boundary conditions: Eigenvalues and eigenvectors
bcs_funcs.bcs_eig = {@bcs_eig};
% bcs_funcs.bcs_eig = bcs_eig_symbolic();

% Boundary conditions: Initial condition
bcs_funcs.bcs_initial = {@bcs_initial};
% bcs_funcs.bcs_initial = bcs_initial_symbolic();

% Boundary conditions: Final condition
bcs_funcs.bcs_final = {@bcs_final};
% bcs_funcs.bcs_final = bcs_final_symbolic();

% Lin gap boundary conditions
bcs_funcs.bcs_lingap = {@bcs_lingap};
% bcs_funcs.bcs_lingap = bcs_lingap_symbolic();

%=========================================================================%
%%                         Initial Continuation                          %%
%=========================================================================%
% We set up the continuation problem by following the equilibrium point x0.
% We then follow the branching point, from which many of the other 
% bifurcations come about.

%-------------------------------------------------------------------------%
%%                    Initial Continuation (Varying A)                   %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.initial_run = 'run01_initial_run';
run_new = run_names.initial_run;

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Initialisation: First Run\n');
fprintf(' Initial continuation from equilibrium point x0\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Continuation parameters : %s\n', 'A');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set upper bound of continuation steps in each direction along solution
prob = coco_set(prob, 'cont', 'PtMX', 50);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);
% Detect and locate Bogdanvos-Takens point
prob = coco_set(prob, 'ep', 'BTP', true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up isol2ep problem
prob = ode_isol2ep(prob, '', funcs.field{:}, x0, pnames, p0);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'gamma'};
prange = {A_range, gamma_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                     Continue From Branching Point                     %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.branching_point = 'run02_branching_point_continuation';
run_new = run_names.branching_point;
% Which run this continuation continues from
run_old = run_names.initial_run;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'BP');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Initialisation: Second Run\n');
fprintf(' Continue bifurcations from the branching point\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set upper bound of continuation steps in each direction along solution
PtMX = 50;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);
% Detect and locate Bogdanvos-Takens point
prob = coco_set(prob, 'ep', 'BTP', true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
% prob = ode_BP2ep(prob, '', run_old, label_old);
prob = ode_ep2ep(prob, '', run_old, label_old);
prob = coco_set(prob, 'cont', 'branch', 'switch');

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'gamma'};
prange = {A_range, gamma_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                         Hopf Bifurcations (H)                         %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.hopf_bifurcations = 'run03_hopf_bifurcation_line_H';
run_new = run_names.hopf_bifurcations;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'HB');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Hopf Bifurcation (H)\n');
fprintf(' Calculate line of Hopf bifurcations\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, 'gamma');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2);
prob = coco_set(prob, 'cont', 'h0', 1e-2);
prob = coco_set(prob, 'cont', 'h_max', 1e-2);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set upper bound of continuation steps in each direction along solution
PtMX = 1000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set number of points
prob = coco_set(prob, 'coll', 'NTST', 700);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
prob = ode_HB2HB(prob, '', run_old, label_old);

%-------------------------%
%     Add COCO Events     %
%-------------------------%
% Add saved point for double limit cycle calculation
prob = coco_add_event(prob, 'DL_PT', 'A', 6.715);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'gamma'};
prange = {A_range, gamma_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                           Saddle-Node (A_S)                           %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.saddle_nodes = 'run04_saddle_node_line_AS';
run_new = run_names.saddle_nodes;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SN');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Saddle Nodes\n');
fprintf(' Calculate line of saddle-nodes\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, gamma');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set upper bound of continuation steps in each direction along solution
PtMX = 20;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from saddle-node
prob = ode_SN2SN(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'gamma'};
prange = {A_range, gamma_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                          Transcritical (A_T)                          %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.transcritical = 'run05_transcritical_points_line_AT';
run_new = run_names.transcritical;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = label_old(3);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Transcritical Points (A_T)\n');
fprintf(' Calculate line of transcritical points A_T\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'gamma, A');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set upper bound of continuation steps in each direction along solution
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Detect and locate neutral saddles
prob = coco_set(prob, 'ep', 'NSA', true);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
prob = ode_ep2ep(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'gamma', 'A'};
prange = {gamma_range, A_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%=========================================================================%
%%                      Compute Double Limit Cycle                       %%
%=========================================================================%
% Compute and follow the double limit cycles (saddle-nodes of periodic
% orbits).

%-------------------------------------------------------------------------%
%%                   Calculate Initial Periodic Orbit                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.limit_cycle.initial_PO = 'run06_initial_periodic_orbit';
run_new = run_names.limit_cycle.initial_PO;
% Which run this continuation continues from
run_old = run_names.hopf_bifurcations;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DL_PT');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Double Limit Cycle: First Run\n');
fprintf(' Compute from Hopf to find saddle-node periodic orbit\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'gamma, A');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
% prob = coco_set(prob, 'cont', 'norm', inf);

% Set NTST mesh
prob = coco_set(prob, 'coll', 'NTST', 25);

% The value of 10 for 'NAdapt' implied that the trajectory discretisation
% is changed adaptively ten times before the solution is accepted.
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set upper bound of continuation steps in each direction along solution
PtMX = 150;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
prob = ode_HB2po(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'gamma', 'A'};
prange = {gamma_range, A_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);


%-------------------------------------------------------------------------%
%%                Continue Saddle Node of Periodic Orbits                %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.limit_cycle.follow_limit_cycle = 'run07_double_limit_cycle';
run_new = run_names.limit_cycle.follow_limit_cycle;
% Which run this continuation continues from
run_old = run_names.limit_cycle.initial_PO;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'SN');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Double Limit Cycle: Second Run\n');
fprintf(' Continue saddle-node of periodic orbits\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, gamma');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 2.5e0);
prob = coco_set(prob, 'cont', 'h0', 2.5e0);
prob = coco_set(prob, 'cont', 'h_max', 2.5e0);

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set upper bound of continuation steps in each direction along solution
PtMX = 5000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set NPR to save every 100 steps
prob = coco_set(prob, 'cont', 'NPR', 100);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
prob = ode_SN2SN(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'gamma'};
prange = {A_range, gamma_range};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%=========================================================================%
%%                Compute Homoclinic Orbits (Approximate)                %%
%=========================================================================%
% We compute a family of homoclinic orbits emanating from a Hopf
% bifurcation. We approximate the homoclinic orbit as a periodic orbit with
% very large period.

%-------------------------------------------------------------------------%
%%        Compute Family of Periodic Orbits from Hopf Bifurcation        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.approx_homo.PO_from_hopf = 'run08_PO_from_hopf_bifurcation';
run_new = run_names.approx_homo.PO_from_hopf;
% Which run this continuation continues from
run_old = run_names.branching_point;

% Label for Hopf bifurcation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'HB');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Approximate Homoclinic: First Run\n');
fprintf(' Continue periodic orbits from a Hopf bifurcation\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, po.period');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detections
prob = coco_set(prob, 'po', 'bifus', 'off');

% Set NAdapt to 1?
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set upper bound of continuation steps in each direction along solution
% 'PtMX', [negative steps, positive steps]
PtMX = 500;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from branching point
prob = ode_HB2po(prob, '', run_old, label_old);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'A', 'po.period'};
prange = {A_range, [0, 1e6]};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                        Large T Periodic Orbit                         %%
%-------------------------------------------------------------------------%
% We approximate a homoclinic orbit by finding a periodic orbit that has a
% large period. The closer you are to an equilibrium point, the longer you
% spend then, and hence you have a super duper large period.

% We obtain an initial solution guess by extending the duration spent near
% an equilibrium point for a periodic orbit found in the previous run.
% The call to coco_xchg_pars constrains the interval length, while
% releasing the second problem parameter.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.approx_homo.large_period_PO = 'run09_reconstruct_large_T_PO';
run_new = run_names.approx_homo.large_period_PO;
% Which run this continuation continues from
run_old = run_names.approx_homo.PO_from_hopf;

% Read solution of previous run with largest period.
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Approximate Homoclinic: Second Run\n');
fprintf(' Find reconstructed high-period periodic orbit approximating a homoclinic\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, homo.po.orb.coll.err_TF, homo.po.period');
fprintf(' =====================================================================\n');

%-------------------%
%     Read Data     %
%-------------------%
% Find minimum point corresponding to equilibrium point, and insert
% large time segment.
data_PO_isol = insert_large_time_segment(run_old, label_old, funcs.field);

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set NTST from previous run
prob = coco_set(prob, 'coll', 'NTST', data_PO_isol.NTST);

% Turn off bifurcation detection
prob = coco_set(prob, 'po', 'bifus', 'off');

% The value of 10 for 'NAdapt' implied that the trajectory discretisation
% is changed adaptively ten times before the solution is accepted.
prob = coco_set(prob, 'cont', 'NAdapt', 10);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue periodic orbit from initial solution
prob = ode_isol2po(prob, 'homo', funcs.field{:}, ...
                   data_PO_isol.tbp_extend, data_PO_isol.xbp_extend, ...
                   data_PO_isol.pnames, data_PO_isol.p);

% Add instance of equilibrium point to find and follow the actual 
% equilibrium point
prob = ode_isol2ep(prob, 'x0', funcs.field{:}, data_PO_isol.x0, data_PO_isol.p);

% Glue parameters
prob = glue_parameters(prob);

%------------------%
%     Run COCO     %
%------------------%
% Assign 'gamma' to the set of active continuation parameters and 'po.period'
% to the set of inactive continuation parameters, thus ensuring that the
% latter is fixed during root finding.
prob = coco_xchg_pars(prob, 'gamma', 'homo.po.period');

% Parameter range
prange = {A_range, [], []};

% Run COCO continuation
coco(prob, run_new, [], 0, {'A', 'homo.po.orb.coll.err_TF', 'homo.po.period'}, prange);

%-------------------------------------------------------------------------%
%%                    Follow Family of "Homoclinics"                     %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.approx_homo.continue_homoclinics = 'run10_approximate_homoclinic';
run_new = run_names.approx_homo.continue_homoclinics;
% Which run this continuation continues from
run_old = run_names.approx_homo.large_period_PO;

% Grab the label for the previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'EP');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Approximate Homoclinic: Third Run\n');
fprintf(' Continue family of periodic orbits approximating homoclinics\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, gamma, homo.po.period');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Set number of continuation steps
PtMX = 10000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Set number of steps to confirm solution
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Turn off bifurcation detections
prob = coco_set(prob, 'po', 'bifus', 'off');

% Set MXCL to false
prob = coco_set(prob, 'coll', 'MXCL', false);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue a periodic orbit from a previous periodic orbit
prob = ode_po2po(prob, 'homo', run_old, label_old);
% Continue from equilibrium point
prob = ode_ep2ep(prob, 'x0', run_old, label_old);

% Glue parameters
prob = glue_parameters(prob);

%------------------%
%     Run COCO     %
%------------------%
% Assign 'gamma' to the set of active continuation parameters and 'po.period'
% to the set of inactive continuation parameters, thus ensuring that the
% latter is fixed during root finding.
prob = coco_xchg_pars(prob, 'gamma', 'homo.po.period');

% Set continuation parameters and range
pcont = {'A', 'gamma', 'homo.po.period'};
prange = {A_range, gamma_range, []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%=========================================================================%
%%                Compute Homoclinic Orbits (Lin's Method)               %%
%=========================================================================%
% We use Lin's method to calculate the homoclinic orbits. This is a better
% approximation than the large period periodic orbits, as above.

%-------------------------------------------------------------------------%
%%               Grow Unstable Manifold of Stationary Point              %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.unstable_manifold = 'run11_unstable_manifold';
run_new = run_names.lins_method.unstable_manifold;
% Which run this continuation continues from
run_old = run_names.approx_homo.continue_homoclinics;

% Grab the label for the previous run solution
label_old = 1;

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: First Run\n");
fprintf(' Continue unstable trajectory segment until we hit Sigma plane\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'seg_u, T1, T2');
fprintf(' =====================================================================\n');

%-------------------%
%     Read Data     %
%-------------------%
% Initial Lin's Method data structure
data_lins = calc_initial_solution_lins(run_old, label_old);

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NTST size
prob = coco_set(prob, 'coll', 'NTST', 25);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set Continuation steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Construct first instance of 'coll' toolbox for unstable manifold
prob = ode_isol2coll(prob, 'unstable', funcs.field{:}, ...
                     data_lins.t0, data_lins.x_init_u, data_lins.pnames, data_lins.p0);
% Construct second instance of 'coll' toolbox for stable manifold
prob = ode_isol2coll(prob, 'stable', funcs.field{:}, ...
                     data_lins.t0, data_lins.x_init_s, data_lins.p0);

% Construct instance of 'ep' tool box to follow stationary point x_neg
% for initial conditions.
prob = ode_isol2ep(prob, 'xneg', funcs.field{:}, data_lins.xneg, data_lins.p0);
% Construct instance of 'ep' tool box to follow stationary point x_pos
% for Lin's segment conditions.
prob = ode_isol2ep(prob, 'xpos', funcs.field{:}, data_lins.xpos, data_lins.p0);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lins, bcs_funcs, false);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'seg_u', 'T1', 'T2'};
prange = {[], [0, 1e3], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                Grow Stable Manifold of Stationary Point               %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.stable_manifold = 'run12_stable_manifold';
run_new = run_names.lins_method.stable_manifold;
% Which run this continuation continues from
run_old = run_names.lins_method.unstable_manifold;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DelU');
label_old = sort(label_old);
label_old = label_old(2);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: Second Run\n");
fprintf(' Continue stable trajectory segment until we hit Sigma plane\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'seg_s, T2, T1');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

% Set Continuation steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue frp,m previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
% Continue frp,m previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Read epsilon and eig data from previous run
data_lingap = read_data_lins(run_old, label_old, false);

% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lingap, bcs_funcs, false);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'seg_s', 'T2', 'T1'};
prange = {[-20, 0], [0, 1e3], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                          Closing the Lin Gap                          %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.close_lingap = 'run13_close_lingap';
run_new = run_names.lins_method.close_lingap;
% Which run this continuation continues from
run_old = run_names.lins_method.stable_manifold;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'DelS');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: Third Run\n");
fprintf(' Close the Lin Gap on the Sigma Plane\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'lingap, T1, T2, theta, A, seg_u');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set Continuation steps
PtMX = 500;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue frp,m previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
% Continue frp,m previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Read epsilon and eig data from previous run
data_lingap = read_data_lins(run_old, label_old, false);

% Calculate Lin gap and vector
data_lingap = calc_lingap_vector(run_old, label_old, data_lingap);

% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lingap, bcs_funcs, true);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'lingap', 'T1', 'T2', 'theta', 'A', 'seg_u'};
% pcont = {'lingap', 'eps1', 'eps2', 'theta', 'gamma', 'seg_u'};
prange = {[0, data_lingap.lingap], [], [], [], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                        Close the Distance eps2                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.close_eps1 = 'run14_close_eps1';
run_new = run_names.lins_method.close_eps1;
% Which run this continuation continues from
run_old = run_names.lins_method.close_lingap;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), 'Lin0');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: Fourth Run\n");
fprintf(' Close epsilon gap until eps1=1e-8\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'eps1, theta, gamma, T1, T2, seg_u');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
% prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set step size
% prob = coco_set(prob, 'cont', 'h_min', 5e-2);
% prob = coco_set(prob, 'cont', 'h0', 5e-2);
% prob = coco_set(prob, 'cont', 'h_max', 2e-1);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set Continuation steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue frp,m previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
% Continue frp,m previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Read epsilon, eig, and lingap data from previous run
data_lingap = read_data_lins(run_old, label_old, true);

% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lingap, bcs_funcs, true);

%------------------------%
%     Add COCO Event     %
%------------------------%
% Add event for when eps1 gets small enough
prob = coco_add_event(prob, 'EPS1', 'eps1', 1e-05);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
% pcont = {'eps1', 'T1', 'theta', 'eps2', 'seg_u', 'gamma'};
pcont = {'eps1', 'theta', 'gamma', 'T1', 'T2', 'seg_u'};
prange = {[1e-8, data_lingap.epsilon(1)], [], [], [], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                        Close the Distance eps2                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.close_eps2 = 'run15_close_eps2';
run_new = run_names.lins_method.close_eps2;
% Which run this continuation continues from
run_old = run_names.lins_method.close_eps1;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), '');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: Fifth Run\n");
fprintf(' Close epsilon gap until eps2=1e-8\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'eps2, theta, gamma, T1, T2, seg_u');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
% prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set step size
% prob = coco_set(prob, 'cont', 'h_min', 5e-2);
% prob = coco_set(prob, 'cont', 'h0', 5e-2);
% prob = coco_set(prob, 'cont', 'h_max', 2e-1);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set Continuation steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue frp,m previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
% Continue frp,m previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Read epsilon, eig, and lingap data from previous run
data_lingap = read_data_lins(run_old, label_old, true);

% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lingap, bcs_funcs, true);

%------------------------%
%     Add COCO Event     %
%------------------------%
% Add event for when eps1 gets small enough
prob = coco_add_event(prob, 'EPS2', 'eps2', 2.0614e-05);

%------------------%
%     Run COCO     %
%------------------%
% Set continuation parameters and range
pcont = {'eps2', 'T2', 'theta', 'eps1', 'seg_u', 'gamma'};
% pcont = {'eps2', 'theta', 'gamma', 'T1', 'T2', 'seg_u'};
prange = {[0, data_lingap.epsilon(2)], [], [], [], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                     Parametrise the Heteroclinic                      %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.lins_method.continue_homoclinics = 'run16_continue_homoclinics';
run_new = run_names.lins_method.continue_homoclinics;
% Which run this continuation continues from
run_old = run_names.lins_method.close_eps2;

% Label for previous run solution
label_old = coco_bd_labs(coco_bd_read(run_old), '');
label_old = max(label_old);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(" Lin's Method: Sixth Run\n");
fprintf(' Continue constrained segments to find parametrisation of homoclinic\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A, gamma, eps1, eps2, theta, seg_u');
fprintf(' =====================================================================\n');

%-------------------------------%
%     Continuation Settings     %
%-------------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set norm
% prob = coco_set(prob, 'cont', 'norm', inf);

% Turn off bifurcation detection
prob = coco_set(prob, 'coll', 'bifus', 'off');

% Turn off MXCL?
prob = coco_set(prob, 'coll', 'MXCL', false);

% Set NAdpat
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% % Set step sizes
% prob = coco_set(prob, 'cont', 'h_min', 1e-1);
% prob = coco_set(prob, 'cont', 'h0', 1e0);
% prob = coco_set(prob, 'cont', 'h_max', 5e0);

% Set Continuation steps
PtMX = 2000;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

% Set frequency of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 50);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Continue from previous solution's segment for the unstable manifold.
prob = ode_coll2coll(prob, 'unstable', run_old, label_old);
% Continue from previous solution's segment for the stable manifold.
prob = ode_coll2coll(prob, 'stable', run_old, label_old);

% Continue frp,m previous 'ep' solution for the x_neg equilibrium point.
prob = ode_ep2ep(prob, 'xneg', run_old, label_old);
% Continue frp,m previous 'ep' solution for the x_pos equilibrium point.
prob = ode_ep2ep(prob, 'xpos', run_old, label_old);

%-----------------------------------%
%     Apply Boundary Conditions     %
%-----------------------------------%
% Read epsilon, eig, and lingap data from previous run
data_lingap = read_data_lins(run_old, label_old, true);

% Glue that shit together, haumi ;)
prob = apply_boundary_conditions_lins(prob, data_lingap, bcs_funcs, true);

%------------------%
%     Run COCO     %
%------------------%
% Calculate Bogdanov Takens Point
gamma_BT = -(1 + (B * (a - 1)) - 2 * sqrt(B * (a - 1))) / (sqrt(B * (a - 1)) * (1 - a - sqrt(B * (a - 1))));

% Set continuation parameters and range
% pcont = {'A', 'gamma', 'T1', 'T2', 'theta', 'seg_u'};
pcont = {'A', 'gamma', 'eps1', 'eps2', 'theta', 'seg_u'};
prange = {A_range, [0, gamma_BT], [], [], [], []};

% Run COCO continuation
coco(prob, run_new, [], 1, pcont, prange);

%=========================================================================%
%%                          SAVE AND PLOT DATA                           %%
%=========================================================================%
%------------------%
%    Save Data     %
%------------------%
% Save data to .mat file
save_fig1_data(run_names, '../data_files/fig1_data.mat');

%----------------------%
%     Plot Figures     %
%----------------------%
% Run plotting scripts
plot_fig1b;
plot_fig1b_inset;

%=========================================================================%
%                               END OF FILE                               %
%=========================================================================%