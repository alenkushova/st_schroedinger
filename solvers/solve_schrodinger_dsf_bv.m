% SOLVE_SCHRODINGER_dsf: Solve a Schrodinger problem with a B-spline
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
%  [geometry, msh, space, u] = solve_schrodinger_dsf (problem_data, method_data)
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

function [geometry, msh, space, u] = solve_schrodinger_dsf_bv (problem_data, method_data)
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
n = numel(nsub);

%trial space-time knots 
[knots, zeta]= kntrefine(geometry.nurbs.knots, nsub-1, trial_degree, trial_regularity);
knots = kntunclamp(knots, trial_degree, trial_regularity, prdc_sides);

[x_knots, x_zeta]= kntrefine(x_geo.nurbs.knots, nsub(1:n-1)-1,...
                             trial_degree(1:n-1), trial_regularity(1:n-1));
x_knots = kntunclamp(x_knots, trial_degree(1:n-1), trial_regularity(1:n-1), prdc_sides);

[t_knots, t_zeta]= kntrefine(t_geo.nurbs.knots,...
                    nsub(n)-1, trial_degree(n), trial_regularity(n));
t_knots = kntunclamp(t_knots, trial_degree(n), trial_regularity(n), prdc_sides);

%test space-time knots
[test_knots, test_zeta]= kntrefine(geometry.nurbs.knots, nsub-1, test_degree, test_regularity);
test_knots = kntunclamp(test_knots, test_degree, test_regularity, prdc_sides);

[test_x_knots, test_x_zeta]= kntrefine(x_geo.nurbs.knots, nsub(1:n-1)-1,...
                             test_degree(1:n-1), test_regularity(1:n-1));
test_x_knots = kntunclamp(test_x_knots, test_degree(1:n-1), test_regularity(1:n-1), prdc_sides);

[test_t_knots, test_t_zeta]= kntrefine(t_geo.nurbs.knots,...
                    nsub(n)-1, test_degree(n), test_regularity(n));
test_t_knots = kntunclamp(test_t_knots, test_degree(n), test_regularity(n), prdc_sides);


% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

[xqn, xqw] = msh_set_quad_nodes (x_zeta, rule(1:n-1));
x_msh   = msh_cartesian (x_zeta, xqn, xqw, x_geo);

[tqn, tqw] = msh_set_quad_nodes (t_zeta, rule(n));
t_msh   = msh_cartesian (t_zeta, tqn, tqw, t_geo);

% Construct msh structures
space   = sp_bspline (knots, trial_degree, msh);
x_space = sp_bspline (x_knots, trial_degree(1:n-1), x_msh);
t_space = sp_bspline (t_knots, trial_degree(n), t_msh);
test_space   = sp_bspline (test_knots, test_degree, msh);
x_test_space = sp_bspline (test_x_knots, test_degree(1:n-1), x_msh);
t_test_space = sp_bspline (test_t_knots, test_degree(n), t_msh);

% Assemble matrices
Wt = op_gradu_v_tp (t_space, t_test_space, t_msh); 
Wt = Wt.';
Mt = op_u_v_tp (t_space, t_test_space, t_msh);
Ms = op_u_v_tp (x_space, x_test_space, x_msh);
Ks = op_laplaceu_v_tp (x_space, x_test_space, x_msh);

% Debug test for OP_LAPLACEU_V_TP:
% diff_op1 = op_geom_exterior ({x_knots}, trial_degree(1:n-1));
% diff_op2 = op_geom_exterior ({x_knots(2:end-1)}, trial_degree(1:n-1)-1);
% D1 = diff_op1{1};
% D2 = diff_op2{1};
% test_Ms = op_u_v_tp (x_test_space, x_test_space, x_msh);
% new_Ks = test_Ms*D2*D1;
% max(max(abs(Ks-new_Ks)))
% This has shown an error of the kind 2.8422e-14. 

% we need to avoid this matrix implementation:
A = gmm*kron(Wt,Ms) - eta*kron(Mt,Ks);
F = op_f_v_tp (space, msh, f); %here is the projection of f.

% Apply Dirichlet bpoundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
[~, boundary_dofs] = sp_drchlt_l2_proj (test_space, msh, h, boundary_sides);
u(drchlt_dofs) = u_drchlt;
trial_int_dofs = setdiff (1:space.ndof, drchlt_dofs);
test_int_dofs  = setdiff (1:test_space.ndof, boundary_dofs);
F(trial_int_dofs) = F(trial_int_dofs) - A(test_int_dofs,drchlt_dofs)*u(drchlt_dofs);
u(trial_int_dofs) = A(test_int_dofs,trial_int_dofs)\F(trial_int_dofs); 
end