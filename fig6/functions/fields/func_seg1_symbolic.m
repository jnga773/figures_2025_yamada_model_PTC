function F_coco_out = func_seg1_symbolic()
  % F_coco_out = func_seg1_symbolic()
  %
  % Creates a CoCo-compatible function encoding for the first
  % segment of the phase-resetting problem.
  %
  % Segment 1 goes from gamma_0 to theta_new.

  % State space dimension
  xdim = 3;

  %---------------%
  %     Input     %
  %---------------%
  % State-space variables
  xvec = sym('x', [xdim, 1]);

  % Adjoint equation variables
  wvec = sym('w', [xdim, 1]);

  % System parameters
  syms gam A B a
  p_sys = [gam; A; B; a];

  % Phase resetting parameters
  syms T k theta_old theta_new
  syms mu_s eta
  syms A_perturb theta_perturb phi_perturb
  p_PR = [T; k; theta_old; theta_new;
          mu_s; eta;
          A_perturb; theta_perturb; phi_perturb];

  % Total vectors
  uvec = [xvec; wvec];
  pvec = [p_sys; p_PR];

  %--------------------------%
  %     Calculate Things     %
  %--------------------------%
  % Vector field
  F_vec = yamada_symbolic_field(xvec, p_sys);

  % Vector equation
  vec_eqn = T * theta_new * F_vec;

  % Calculate tranpose of Jacobian at point xvec
  J_T = transpose(jacobian(F_vec, xvec));

  % Adjoint equation
  adj_eqn = -T * theta_new * J_T * wvec;

  % Total equation
  F_seg = [vec_eqn; adj_eqn];

  %-----------------%
  %     SymCOCO     %
  %-----------------%
  % Filename for output functions
  filename_out = './functions/symcoco/F_seg1';

  % COCO Function encoding
  F_coco = sco_sym2funcs(F_seg, {uvec, pvec}, {'x', 'p'}, 'filename', filename_out);

  % List of functions
  func_list = {F_coco(''), ...
               F_coco('x'), F_coco('p'), ...
               F_coco({'x', 'x'}), F_coco({'x', 'p'}), F_coco({'p', 'p'})};

  %----------------%
  %     Output     %
  %----------------%
  F_coco_out = func_list;

end