% This function computes the mass matrix associated to a spline space in
% order to represent the square l2-norm associated to the inner product
%
% Call: M = massmatrix(geometry, msh, space)
% 
function M = massmatrix(geo_name, trial_degree, test_degree, nel)
geometry = geo_load('geo_square.txt');

% Choice of discretization parameters. 
trial_reg = degree-1;                                                      % Regularity of the splines in trial space 
test_reg  = degree-1;                                                      % Regularity of the splines in test space
nsub      = [nel nel];                                                     % Number of subdivisions ( " )
nquad     = max (trial_degree, test_degree) + 1;                           % Points for the Gaussian quadrature rule ( " )

% knots in space-time:
[trial_knots, trial_zeta]= kntrefine (geometry.nurbs.knots, nsub-1,...
                                      trial_degree, trial_reg);
trial_knots = kntunclamp(trial_knots, trial_degree, trial_reg, []);

% Construct msh structures.
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (trial_zeta, rule);
msh      = msh_cartesian (trial_zeta, qn, qw, geometry);

% Construct the space structures.
trial_space = sp_bspline (trial_knots, degree, msh);
test_space  = sp_bspline (test_knots, degree, msh);

M  = op_u_v_tp (space, space, msh);
end