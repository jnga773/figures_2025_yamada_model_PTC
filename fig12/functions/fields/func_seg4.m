function y_out = func_seg4(x_in, p_in)
  % y_out = func_seg1(u_in, p_in)
  %
  % COCO 'ode' toolbox encoding for the vector field corresponding to
  % segment 4 of the phase resetting curve.
  %
  % Segment 4 goes from theta_old to theta_new.
  %
  % Parameters
  % ----------
  % x_in : array, float
  %     State vector for the periodic orbit (x) and perpendicular
  %     vector (w).
  % p_in : array, float
  %     Array of parameter values
  %
  % Returns
  % -------
  % y_out : array, float
  %     Array of the vector field of the periodic orbit segment
  %     and the corresponding adjoint equation for the perpendicular
  %     vector.

  % Original vector field dimensions (CHANGE THESE)
  xdim = 3;
  pdim = 4;
  % Original vector field function
  field      = @yamada;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Array of state-space variables
  x_vec        = x_in(1:xdim, :);

  % System parameters
  p_system     = p_in(1:pdim, :);

  % Phase resetting parameters
  % Integer for period
  k             = p_in(pdim+1, :);
  % Phase where perturbation starts
  % theta_old     = p_in(pdim+2, :);
  % Phase where segment comes back to \Gamma
  % theta_new     = p_in(pdim+3, :);
  % Stable Floquet eigenvalue
  % mu_s          = p_in(pdim+4, :);
  % Distance from pertured segment to \Gamma
  % eta           = p_in(pdim+5, :);
  % Size of perturbation
  % A_perturb     = p_in(pdim+6, :);
  % Angle of perturbation
  % theta_perturb = p_in(pdim+7, :);
  % Azimuthal angle of perturbation
  % phi_perturb   = p_in(pdim+8, :);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%  
  % Calculate vector field
  vec_field = field(x_vec, p_system);
  
  % Save to array
  vec_eqn = k .* vec_field;

  %----------------%
  %     Output     %
  %----------------%
  % Vector field
  y_out(1:xdim, :) = vec_eqn(:, :);

end
