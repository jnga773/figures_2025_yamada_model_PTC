function y_out = VAR(u_in, p_in)
  % y_out = VAR(u_in, p_in)
  %
  % COCO 'ode' toolbox encoding for the adjoint problem, solving for
  % the rotated periodi orbit and the stable Floquet bundle.
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
  field_DFDX = @yamada_DFDX;

  %---------------%
  %     Input     %
  %---------------%
  % State vector
  x          = u_in(1 : xdim, :);

  % Adjoint perpendicular vector (I think)
  w          = u_in(xdim+1 : 2*xdim, :);

  % Parameters
  parameters = p_in(1:end, :);
  
  % System parameters
  p_system = parameters(1:pdim, :);

  % Floquet eigenvalue
  mu_s  = parameters(pdim+1, :);
  % Norm w-vector
  wnorm = parameters(pdim+2, :);

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field
  vec_field = field(x, p_system);

  % Jacobian
  J = field_DFDX(x, p_system);

  % Vector field equations
  vec_eqn = vec_field;

  % Length of u_in
  ll = length(x(1, :));

  % Empty array for adjoint equations
  adj_eqn = zeros(xdim, ll);

  % The adjoint equation
  % d/dt w = -J^T(\gamma) * w
  for i = 1 : ll
    % Transpose of the Jacobian
    J_transpose = transpose(J(:, :, i));

    % Adjoint equation
    adj_eqn(:, i) = -J_transpose * w(:, i);

  end

  %----------------%
  %     Output     %
  %----------------%  
  % Vector field
  y_out(1:xdim, :) = vec_eqn(:, :);


  % Adjoint equation
  y_out(xdim+1:2*xdim, :) = adj_eqn(:, :);

end
