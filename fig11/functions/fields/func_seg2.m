function y_out = func_seg2(x_in, p_in)
  % y_out = func_seg2(x_in, p_in)
  %
  % Creates a CoCo-compatible function encoding for the second
  % segment of the phase-resetting problem.
  %
  % Segment 2 goes from \gamma_{\vartheta} to \gamma_{0}.
  %
  % Parameters
  % ----------
  % x_in : array, double
  %     State vector for the periodic orbit (x) and perpendicular
  %     vector (w).
  % p_in : array, double
  %     Array of parameter values
  %
  % Returns
  % -------
  % y_out : array, double
  %     Array of the vector field of the periodic orbit segment
  %     and the corresponding adjoint equation for the perpendicular
  %     vector.

  %============================================================================%
  %                          CHANGE THESE PARAMETERS                           %
  %============================================================================%
  % Original vector field state-space dimension
  xdim       = 3;
  % Original vector field parameter-space dimension
  pdim       = 4;
  % Original vector field function
  field      = @yamada;

  %============================================================================%
  %                                    INPUT                                   %
  %============================================================================%
  %-------------------------------%
  %     State-Space Variables     %
  %-------------------------------%
  % State space variables
  x_vec = x_in(1:xdim, :);
  % Perpendicular vectors
  w_vec = x_in(xdim+1:2*xdim, :);
  
  %--------------------%
  %     Parameters     %
  %--------------------%
  % System parameters
  p_sys = p_in(1:pdim, :);
  % Phase along \Gamma
  theta = p_in(pdim+1, :);

  %============================================================================%
  %                           VECTOR FIELD ENCODING                            %
  %============================================================================%
  %----------------------%
  %     Vector Field     %
  %----------------------%
  % Calculate vector field
  vec_field = field(x_vec, p_sys);
  
  % Save to array
  vec_eqn = (1.0 - theta) .* vec_field;

  %============================================================================%
  %                                   OUTPUT                                   %
  %============================================================================%
  % Vector field
  y_out(1:xdim, :) = vec_eqn(:, :);

end
