% SOLVE_SCHRODINGER_ST_New: Solve a Schrodinger problem with a B-spline
%                               discretization in a space-time isoparametric
%                               approach . 
%
% The function solves the evolution problem
%
%     i d_t (u) - d_x( d_x (u)) = f    in Omega   = (0,1)^d x [0,T] 
%                             u = g    on Omega_0 = (0,1)^d x {0}
%                             u = h    on Gamma_D
%
% using the discrete strong formulation (well posed problem).
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_schrodinger (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition 
%    - prdc_sides:   sides with Periodic boundary condition (may be empty)
%    - f:            source term
%    - g:            function for Neumann boundary condition (if needed)
%    - h:            function for Dirichlet boundary condition (also initial)
%    - T:            final time T = 10 by default specified in Geo file
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions
%    - regularity: continuity of the spline functions
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - solver:     Whether 'FD' for Fast Diagonalization or 'M' Matlab backslash 
%                  (Defoult Matlab backslash)
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete space (see sp_scalar)
%  u:        the computed degrees of freedom
%
% See also EX_SCHRODINGER_A or EX_SCHRODINGER_B for examples.
%

function [geometry, msh, space, u] = solve_schrodinger(problem_data, method_data)
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end
% Construct geometry structures
geometry = geo_load(xt_geo_name);
x_geo    = geo_load(x_geo_name);
t_geo    = geo_load(t_geo_name);

%trial space-time knots 
[knots, zeta]= kntrefine(geometry.nurbs.knots, nsub-1, trial_degree, trial_regularity);
knots = kntunclamp(knots, trial_degree, trial_regularity, prdc_sides);

%test space-time knots
[test_knots, test_zeta]= kntrefine(geometry.nurbs.knots, nsub-1, test_degree, test_regularity);
test_knots = kntunclamp(test_knots, test_degree, test_regularity, prdc_sides);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

[xqn, xqw] = msh_set_quad_nodes (zeta(1), rule(1));
x_msh   = msh_cartesian (zeta(1), xqn, xqw, x_geo);

[tqn, tqw] = msh_set_quad_nodes (zeta(2), rule(2));
t_msh   = msh_cartesian (zeta(2), tqn, tqw, t_geo);

% Construct msh structures
space        = sp_bspline (knots, trial_degree, msh);
x_space      = sp_bspline (knots(1), trial_degree(1),x_msh);
t_space      = sp_bspline (knots(2), trial_degree(2),t_msh);
test_space   = sp_bspline (test_knots, test_degree, msh);
x_test_space = sp_bspline (test_knots(1), test_degree(1),x_msh);
t_test_space = sp_bspline (test_knots(2), test_degree(2),t_msh);

% Assemble matrices
Wt = op_gradu_v_tp (t_space, t_test_space, t_msh); %Controllo come è fatto!
Wt = Wt';
Mt = op_u_v_tp (t_space, t_test_space, t_msh);
Ms = op_u_v_tp (x_space, x_test_space, x_msh);
Ks = op_laplaceu_v_tp (x_space, x_test_space, x_msh);

% we need to avoid this matrix implementation:
A = gmm*kron(Wt,Ms) - eta*kron(Mt,Ks);
F = op_f_v_tp (space, msh, f); %here is the projection of f.

% Apply Dirichlet bpoundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
F(int_dofs) = F(int_dofs) - A(:,drchlt_dofs)*u(drchlt_dofs);
u(int_dofs) = A(:,int_dofs)\F(int_dofs);
end