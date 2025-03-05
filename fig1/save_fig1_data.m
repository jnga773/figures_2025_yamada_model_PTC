function save_fig1_data(run_names_in, filename_in)
  % save_fig1_data(run_names_in, filename_in)
  %
  % Reads the bifurcation data from all of the runs in [run_names_in], 
  % and saves the data to filename_in.

  %--------------------------------------%
  % Read Run Name String Identifiers     %
  %--------------------------------------%
  H_run        = run_names_in.hopf_bifurcations;
  S_run        = run_names_in.saddle_nodes;
  T_run        = run_names_in.transcritical;
  L_run_approx = run_names_in.approx_homo.continue_homoclinics;
  L_run_lins   = run_names_in.lins_method.continue_homoclinics;
  D_run        = run_names_in.limit_cycle.follow_limit_cycle;
    
  %-------------------------------------%
  %     Read Bifurcation Data Files     %
  %-------------------------------------%
  % Read COCO data matrices
  bd_H        = coco_bd_read(H_run);
  bd_S        = coco_bd_read(S_run);
  bd_T        = coco_bd_read(T_run);
  bd_L_approx = coco_bd_read(L_run_approx);
  bd_L_lins   = coco_bd_read(L_run_lins)
  bd_D        = coco_bd_read(D_run);

  %----------------------------------------%
  %     Read Data: Bifurcation Diagram     %
  %----------------------------------------%
  % Hopf bifurcation line (H)
  A_run8 = coco_bd_col(bd_H, 'A');
  gamma_run8 = coco_bd_col(bd_H, 'gamma');

  % Find the minimum to split into H line and NSA line
  [~, idx] = min(A_run8);

  % Neutral saddle bifurcations
  A_NSA     = A_run8(1:idx)
  gamma_NSA = gamma_run8(1:idx);

  % Hopf bifurcations
  A_H     = A_run8(idx+1:end)
  gamma_H = gamma_run8(idx+1:end);

  % Saddle-Node bifurcation line (A_S)
  A_SN     = coco_bd_col(bd_S, 'A');
  gamma_SN = coco_bd_col(bd_S, 'gamma');

  % Transcritical bifurcation line (A_T)
  A_T     = coco_bd_col(bd_T, 'A');
  gamma_T = coco_bd_col(bd_T, 'gamma');

  % Approximate homoclinic line
  A_L_approx     = coco_bd_col(bd_L_approx, 'A');
  gamma_L_approx = coco_bd_col(bd_L_approx, 'gamma');

  % Lin's method homoclinic line
  A_L_lins     = coco_bd_col(bd_L_lins, 'A');
  gamma_L_lins = coco_bd_col(bd_L_lins, 'gamma');

  % Approximate double limit cycle line
  A_D     = coco_bd_col(bd_D, 'A');
  gamma_D = coco_bd_col(bd_D, 'gamma');

  %----------------%
  %     Output     %
  %----------------%
  % Hopf
  data_out.A_H            = A_H;
  data_out.gamma_H        = gamma_H;

  % Neutral-saddle
  data_out.A_NSA          = A_NSA;
  data_out.gamma_NSA      = gamma_NSA;

  % Saddle-node
  data_out.A_S            = A_S;
  data_out.gamma_S        = gamma_S;

  % Transcritical
  data_out.A_T            = A_T;
  data_out.gamma_T        = gamma_T;

  % Double-limit cycle
  data_out.A_D            = A_D;
  data_out.gamma_D        = gamma_D;

  % Homoclinic: Approximate
  data_out.A_L_approx     = A_L_approx;
  data_out.gamma_L_approx = gamma_L_approx;

  % Homoclinic: Lins
  data_out.A_L_lins       = A_L_lins;
  data_out.gamma_L_lins   = gamma_L_lins;

  %-------------------%
  %     Save Data     %
  %-------------------%
  save(filename_in, '-struct', 'data_out');

end
