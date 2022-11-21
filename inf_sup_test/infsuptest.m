% INF_SUP_TEST for Schroedinger space-time discrete weak formulaiton
%  p   =  degree of B-splines
%  nel =  number of subdivisions in space direction!
% 
% This tests if the infsup discrete condition is achieved.
% 
% we check inf_v sup_w (A*v,w)/ |v||w| 
% A is the matrix associated to the bilinear form.
% Mv is the matrix representing the norm on trial space.
% Mw is the matrix representing the norm on test space.
% 
% It also can compute the eigenfunction associated to the smallestabs
% eigenvalue, in order to visualize its behaviour. 
function [mu, D, geo, msh, space, w] = infsuptest (p, nel)
% Boundaries.
drchlt_sides   = [1 2 3 ];% Dirichlet.
h   = @(x, t, ind) 0*x;
% Set the discretization
[geo, msh, mshx, msht, space, x_space, t_space] = set_schrodinger_dwf (p,nel);

% Assembly the matrices.
A = op_schroedinger_1st_order (mshx, msht, x_space, t_space);
Mw  = op_u_v_tp (space, space, msh);
Mv = op_schroedinger_graph_norm (msh, mshx, msht, space, x_space, t_space);

% Extract internal dofs.
[~, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
A  = A  (int_dofs, int_dofs);
Mw = Mw (int_dofs, int_dofs);
Mv = Mv (int_dofs, int_dofs);

% Compute the inf-sup constants:
%[nu, C] = infsup_coercivity(A, Mv, Mw);
[mu, w, D] = infsup_suriettivity(A, Mv, Mw);
mu = sqrt(real(mu));
end