function [data_in, y_out] = bcs_eig(prob_in, data_in, u_in)
  % [data_in, y_out] = bcs_eig(prob_in, data_in, u_in)
  % 
  % COCO compatible encoding for the boundary conditions of the eigenvalues and
  % eigenvectors of the monodromy matrix. Ensures they are eigenvectors and
  % values of the monodromy matrix, and ensures that the eigenvector is
  % normalised.
  %
  % Input
  % ----------
  % prob_in : COCO problem structure
  %     Continuation problem structure.
  % data_in : structure
  %     Problem data structure contain with function data.
  % u_in : array (floats?)
  %     Total u-vector of the continuation problem. This function
  %     only utilises the following (as imposed by coco_add_func):
  %          * u_in(1:3)  - Equilibrium point
  %          * u_in(4:7)  - System parameters
  %          * u_in(8:10) - The eigenvector
  %          * u_in(11)   - The eigenvalue
  %
  % Output
  % ----------
  % y_out : array of vectors
  %     An array containing to the two boundary conditions.
  % data_in : structure
  %     Not actually output here but you need to have it for COCO.

  % State space dimension
  xdim = data_in.xdim;
  % Parameter space dimension
  pdim = data_in.pdim;
  % Jacobian
  field_DFDX = data_in.dfdxhan;
  % field_DFDX = yamada_DFDX;

  %--------------------------%
  %     Input Parameters     %
  %--------------------------%
  % Equilibrium point
  x0_ss      = u_in(1 : xdim);
  % System parameters
  parameters = u_in(xdim+1 : xdim+pdim);
  % Eigenvector
  eig_vec    = u_in(xdim+pdim+1 : xdim+pdim+xdim);
  % Eigenvalue
  eig_val    = u_in(end);

  %----------------------------%
  %     Calculate Jacobian     %
  %----------------------------%
  % Calculate Jacobian matrix
  % J = yamada_DFDX(x0_ss, parameters);
  J = field_DFDX(x0_ss, parameters);

  %---------------------------------------%
  %     Calculate Boundary Conditions     %
  %---------------------------------------%
  % Eigenvalue equations
  eig_eqn = (J * eig_vec) - (eig_val * eig_vec);

  % Unit vector equations
  vec_eqn = (eig_vec' * eig_vec) - 1;

  %----------------%
  %     Output     %
  %----------------%
  % Boundary conditions
  y_out = [eig_eqn;
           vec_eqn];

end