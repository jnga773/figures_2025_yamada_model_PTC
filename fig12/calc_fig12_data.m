%=========================================================================%
%%               Compute Floquet Bundle at Zero Phase Point              %%
%=========================================================================%
% We now add the adjoint function and Floquet boundary conditions to
% compute the adjoint (left or right idk) eigenvectors and eigenvalues.
% This will give us the perpendicular vector to the tangent of the periodic
% orbit. However, this will only be for the eigenvector corresponding to
% the eigenvalue \mu = 1.

%-------------------------------------------------------------------------%
%%                     Compute Stable Eigenvalue 1.0                     %%
%-------------------------------------------------------------------------%
% Starting from an initial zero vector, we continue in mu until the stable
% eigenvalue is 1.0

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.VAR_mu = 'run03_VAR_mu';
run_new = run_names.VAR_mu;
% Which run this continuation continues from
run_old = run_names.initial_PO_COLL;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'PO_PT');

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Floquet Bundle: First Run\n');
fprintf(' Calculate stable Floquet bundle eigenvalue\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'mu_s, w_norm');
fprintf(' =====================================================================\n');

%--------------------------%
%     Calculate Things     %
%--------------------------%
data_adjoint = calc_initial_solution_VAR(run_old, label_old);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 1e-2, 'h0', 1e-2, 'h_max', 1e-2);

% Set PtMX
PtMX = 100;
prob = coco_set(prob, 'cont', 'PtMX', PtMX);

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
coco(prob, run_new, [], 1, {'mu_s', 'w_norm'} , {[0.9, 1.1], []});

%-------------------------------------------------------------------------%
%%                  Grow Orthogonal Stable Eigenvector                   %%
%-------------------------------------------------------------------------%
% Having found the solution (branching point 'BP') corresponding to
% \mu = 1, we can continue in the norm of the vector w (w_norm), until the
% norm is equal to zero. Then we will have the correct perpendicular
% vector.

%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.VAR_wnorm = 'run04_VAR_wnorm';
run_new = run_names.VAR_wnorm;
% Which run this continuation continues from
run_old = run_names.VAR_mu;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'BP');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Floquet Bundle: Second Run\n');
fprintf(' Grow norm of stable Floquet bundle vector\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'mu_s, w_norm');
fprintf(' =====================================================================\n');

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set number of PtMX steps
PtMX = 200;
prob = coco_set(prob, 'cont', 'PtMX', [0, PtMX]);

% Set number of saved solutions
prob = coco_set(prob, 'cont', 'NPR', 25);

% Continue coll from previous branching point
% prob = ode_BP2coll(prob, 'adjoint', run_old, label_old);
prob = ode_coll2coll(prob, 'adjoint', run_old, label_old);
prob = coco_set(prob, 'cont', 'branch', 'switch');

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
coco(prob, run_new, [], 1, {'mu_s', 'w_norm'}, {[], [-1e-4, 1.1]});

%=========================================================================%
%%                   CALCULATE PHASE RESET SOLUTIONS                     %%
%=========================================================================%
% We compute a set of phase reset problems. We first continue in the
% parameter 'A_perturb', that is, the amplitude of the perturbation.
% After this, we then compute some phase transition curves (PTCs), by
% continuing in 'theta_old' and 'theta_new'.

%-------------------------------------------------------------------------%
%%                    Increase Perturbation Amplitude                    %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.PR_increase_perturbation = 'run05_PR_increase_perturbation';
run_new = run_names.PR_increase_perturbation;
% Which run this PR_increase_perturbation continues from
run_old = run_names.VAR_wnorm;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'NORM1');
label_old = label_old(1);

%--------------------------%
%     Print to Console     %
%--------------------------%
fprintf(' =====================================================================\n');
fprintf(' Directional Transition Curve: First Run\n');
fprintf(' Increase perturbation amplitude\n');
fprintf(' ---------------------------------------------------------------------\n');
fprintf(' This run name           : %s\n', run_new);
fprintf(' Previous run name       : %s\n', run_old);
fprintf(' Previous solution label : %d\n', label_old);
fprintf(' Continuation parameters : %s\n', 'A_perturb, theta_new, eta, mu_s');
fprintf(' =====================================================================\n');

%-------------------%
%     Read Data     %
%-------------------%
% Set periodicity
k = 30;

% Set perturbation direction to be d = (1, 0, 1) / sqrt(2)
% theta_perturb = 0;
theta_perturb = 0.0;

% Set initial conditions from previous solutions
data_PR = calc_initial_solution_PR(run_old, label_old, k, theta_perturb);

%----------------------------%
%     Setup Continuation     %
%----------------------------%
% Set up the COCO problem
prob = coco_prob();

% Set step sizes
prob = coco_set(prob, 'cont', 'h_min', 5e-5);
prob = coco_set(prob, 'cont', 'h0', 1e-2);
prob = coco_set(prob, 'cont', 'h_max', 1e0);

% Set adaptive mesh
prob = coco_set(prob, 'cont', 'NAdapt', 10);

% Set number of steps
prob = coco_set(prob, 'cont', 'PtMX', 5000);

% Set number of stored solutions
prob = coco_set(prob, 'cont', 'NPR', 100);

% Turn off MXCL
prob = coco_set(prob, 'coll', 'MXCL', 'off');

% Set norm to int
prob = coco_set(prob, 'cont', 'norm', inf);

% % Set MaxRes and al_max
% prob = coco_set(prob, 'cont', 'MaxRes', 10);
% prob = coco_set(prob, 'cont', 'al_max', 25);

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
                     data_PR.t_seg1, data_PR.x_seg1, data_PR.pnames, data_PR.p0);

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
% List of perturbation amplitudes to save solutions for
SP_parameter = 'A_perturb';
SP_values    = [0.1, 0.724236, 25.0];

% Save solution at phase along \Gamma where there WILL BE an intersection
% with the stable manifold of q.
prob = coco_add_event(prob, 'SP', SP_parameter, SP_values);

%------------------%
%     Run COCO     %
%------------------%

% Set continuation parameters and parameter range
pcont = {'A_perturb', 'theta_new', ...
         'eta', 'mu_s'};
prange = {[0.0, max(SP_values)+.1], [], ...
          [-1e-4, 1e-2], [0.99, 1.01]};

% Run COCO
coco(prob, run_new, [], 1, pcont, prange);

%-------------------------------------------------------------------------%
%%                      Change Perturbation Angle                        %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.PR_change_angle = 'run06_PR_change_angle';
run_new = run_names.PR_change_angle;
% Which run this continuation continues from
run_old = run_names.PR_increase_perturbation;

% Continuation point
label_old = coco_bd_labs(coco_bd_read(run_old), 'SP');

%---------------------------------%
%     Cycle through SP labels     %
%---------------------------------%
% Set number of threads
% M = 3;
% parfor (run = 1 : length(label_old), M)
for run = 1 : length(label_old)
  % Label for this run
  this_run_label = label_old(run);

  % Data directory for this run
  this_run_name = {run_new; sprintf('A_perturb_%02d', run)};

  %--------------------------%
  %     Print to Console     %
  %--------------------------%
  fprintf(' =====================================================================\n');
  fprintf(' Directional Transition Curve: Second Run\n');
  fprintf(' Change perturbation angles for two approaches\n');
  fprintf(' ---------------------------------------------------------------------\n');
  fprintf(' This run name           : {%s, %s}\n', this_run_name{1}, this_run_name{2});
  fprintf(' Previous run name       : %s\n', run_old);
  fprintf(' Previous solution label : %d\n', this_run_label);
  fprintf(' Continuation parameters : %s\n', 'theta_perturb, theta_new, eta, mu_s');
  fprintf(' =====================================================================\n');

  %------------------%
  %     Run COCO     %
  %------------------%
  % Saved points
  SP_parameter = 'theta_perturb';
  SP_values = [0.0, 0.4];

  if run == 3
    SP_values(2) = 0.28;
  end

  % Continuation parameters
  pcont = {'theta_perturb', 'theta_new', 'eta', 'mu_s', 'A_perturb'};
  % Parameter range for continuation
  prange = {[-1e-3, max(SP_values)+0.01], [], [-1e-4, 1e-2], [0.99, 1.01], []};

  % Run continuation
  run_PTC_continuation(this_run_name, run_old, this_run_label, data_PR, bcs_funcs, ...
                       pcont, prange, ...
                       SP_parameter=SP_parameter, SP_values=SP_values, ...
                       h_min=1e-3, h0=1e-2, h_max=1e1, ...
                       PtMX=[0, 500], NPR=20, NAdapt=20);

end

%-------------------------------------------------------------------------%
%%                     Move theta_old Along \Gamma                       %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.PR_move_theta_old = 'run07_PR_move_theta_old';
run_new = run_names.PR_move_theta_old;
% Which run this continuation continues from
run_old = run_names.PR_change_angle;

% Folder name
dirs_A = sprintf('./data/%s/', run_old);
% List all directories
dirs_A = dir(dirs_A);
% Remove ./ and ../
dirs_A = dirs_A(~ismember({dirs_A.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dirs_A = {dirs_A.name};

%-------------------------------------%
%     Cycle Through Previous Runs     %
%-------------------------------------%
% Set number of threads
N_threads = 3;
% N_threads = length(dirs_A);
parfor(idx = 1 : length(dirs_A), N_threads)
  % Set run string for this run
  sub_run_name = {run_old, dirs_A{idx}};

  % Read bd_data
  bd_read = coco_bd_read(sub_run_name);

  % Read previous labels
  label_old = coco_bd_labs(bd_read, 'SP');

  % Check if there are 2 'SP' solutions or not
  if isscalar(label_old)
    % Append the first EP solution
    EPs = coco_bd_labs(bd_read, 'EP');
    label_old = [min(EPs), label_old];
  end

  % Cycle through labels and move theta_old
  for run = 1 : length(label_old)
    % Label for this run
    this_run_label = label_old(run);

    % Data directory for this run
    this_run_name = {run_new; dirs_A{idx}; sprintf('theta_perturb_%02d', run)};

    %--------------------------%
    %     Print to Console     %
    %--------------------------%
    fprintf(' =====================================================================\n');
    fprintf(' Directional Transition Curve: Third Run\n');
    fprintf(' Move along periodic orbit\n');
    fprintf(' ---------------------------------------------------------------------\n');
    fprintf(' This run name           : {%s, %s, %s}\n', this_run_name{1}, this_run_name{2}, this_run_name{3});
    fprintf(' Previous run name       : {%s, %s}\n', sub_run_name{1}, sub_run_name{2});
    fprintf(' Previous solution label : %d\n', this_run_label);
    fprintf(' Continuation parameters : %s\n', 'theta_old, theta_new, eta, mu_s');
    fprintf(' =====================================================================\n');

    %------------------%
    %     Run COCO     %
    %------------------%
    % Saved points
    SP_parameter = 'theta_old';
    SP_values = [0.339386, 1.339386];

    % Continuation parameters
    pcont = {'theta_old', 'theta_new', 'eta', 'mu_s', 'A_perturb'};
    % Parameter range for continuation
    prange = {[min(SP_values), max(SP_values)], [], [-1e-4, 1e-2], [0.99, 1.01], []};

    % Run continuation
    run_PR_continuation(this_run_name, sub_run_name, this_run_label, data_PR, bcs_funcs, ...
                        pcont, prange, ...
                        SP_parameter=SP_parameter, SP_values=SP_values, ...
                        h_min=1e-3, h0=1e-1, h_max=1e1, ...
                        PtMX=500, NPR=50, NAdapt=20);
  end
end

%-------------------------------------------------------------------------%
%%                       Calculate DTCs (Multiple)                       %%
%-------------------------------------------------------------------------%
%------------------%
%     Run Name     %
%------------------%
% Current run name
run_names.PR_DTC_scan = 'run08_PR_DTC_scan';
run_new = run_names.PR_DTC_scan;
% Which run this continuation continues from
run_old = run_names.PR_move_theta_old;

% Folder name
dirs_A = sprintf('./data/%s/', run_old);
% List all directories
dirs_A = dir(dirs_A);
% Remove ./ and ../
dirs_A = dirs_A(~ismember({dirs_A.name}, {'.', '..', '.DS_Store'}));
% Sub folder names
dirs_A = {dirs_A.name};

%-------------------------------------%
%     Cycle Through Previous Runs     %
%-------------------------------------%
% for idx_A = 1 : length(dirs_A)
parfor (idx_A = 1 : length(dirs_A), 3)
  % Read folders inside dirs_A
  dirs_theta = dir(sprintf('./data/%s/%s/', run_old, dirs_A{idx_A}));
  dirs_theta = dirs_theta(~ismember({dirs_theta.name}, {'.', '..', '.DS_Store'}));
  dirs_theta = {dirs_theta.name};

  % Cycle through angle runs
  for idx_theta = 1 : length(dirs_theta)
    % Set run string for this run
    sub_run_name = {run_old, dirs_A{idx_A}, dirs_theta{idx_theta}};

    % Read bd_data
    bd_read = coco_bd_read(sub_run_name);

    % Read previous labels
    label_old = coco_bd_labs(bd_read, 'SP');

    % Check if label_old is not empty and, if not, run the continuation
    if ~isempty(label_old)
      for run = 1 : length(label_old)
        % Label for this run
        this_run_label = label_old(run);

        % Data directory for this run
        this_run_name = {run_new; dirs_A{idx_A}; dirs_theta{idx_theta}; sprintf('theta_old_%02d', run)};

        %--------------------------%
        %     Print to Console     %
        %--------------------------%
        fprintf(' =====================================================================\n');
        fprintf(' Directional Transition Curve: Third Run\n');
        fprintf(' Move along periodic orbit\n');
        fprintf(' ---------------------------------------------------------------------\n');
        fprintf(' This run name           : {%s, %s, %s, %s}\n', this_run_name{1}, this_run_name{2}, this_run_name{3}, this_run_name{4});
        fprintf(' Previous run name       : {%s, %s, %s}\n', sub_run_name{1}, sub_run_name{2}, sub_run_name{3});
        fprintf(' Previous solution label : %d\n', this_run_label);
        fprintf(' Continuation parameters : %s\n', 'theta_perturb, theta_new, eta, mu_s');
        fprintf(' =====================================================================\n');

        %------------------%
        %     Run COCO     %
        %------------------%
        % Continuation parameters
        pcont = {'theta_perturb', 'theta_new', 'eta', 'mu_s', 'A_perturb'};
        % Parameter range for continuation
        prange = {[0.0, 2.0], [], [-1e-4, 1e-2], [0.99, 1.01], []};

        % Run continuation
        run_PTC_continuation(this_run_name, sub_run_name, this_run_label, data_PR, bcs_funcs, ...
                            pcont, prange, ...
                            h_min=1e-4, h0=5e-2, h_max=1e1, ...
                            PtMX=1000, NPR=50, NAdapt=10);
      end
    end
  end
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