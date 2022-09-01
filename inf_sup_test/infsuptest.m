% INF_SUP_TEST for Schroedinger space-time discrete weak formulaiton
%  p   =  degree of B-splines
%  nel =  number of subdivisions in space direction!
% 
function [mu, D] = infsuptest (p, nel)
% Geometries.
T = 2;
xt_geo_name = 'geo_rectangle.txt';
x_geo_name  = nrbline ([0 0], [1 0]);
t_geo_name  = nrbline ([0 0], [T 0]);

geometry = geo_load(xt_geo_name);
x_geo    = geo_load(x_geo_name);
t_geo    = geo_load(t_geo_name);

% Boundaries.
nmnn_sides     = [];      % Neumann.
drchlt_sides   = [1 2 3 ];% Dirichlet.
x_drchlt_sides = [1 2 ];  % Dirichlet.
prdc_sides     = [];      % Periodic.

% Parameters + exact solution.
omg = 0.2;
a0  = (2/(omg^2))^(1/4);
gmm = 1i;
eta = 1;
uex = @(x, t) (a0*exp(-1i*(x.^2 + t.^2)/(omg^2)));
h   = @(x, t, ind) uex(x,t);

% Choice of discretization parameters. 
degree     = [p p]; % Degree of the splines (last is time dir)
regularity = degree-1; % Regularity of the splines ( " )
nsub       = [nel T*nel]; % Number of subdivisions ( " )
nquad      = degree+1; % Points for the Gaussian quadrature rule ( " )

rdim = numel(degree);

[knots, zeta]= kntrefine(geometry.nurbs.knots, nsub-1, degree, regularity);
knots = kntunclamp(knots, degree, regularity, prdc_sides);

[x_knots, x_zeta]= kntrefine(x_geo.nurbs.knots,...
                    nsub(1:rdim-1)-1, degree(1:rdim-1), regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, degree(1:rdim-1), regularity(1:rdim-1), prdc_sides);

[t_knots, t_zeta]= kntrefine(t_geo.nurbs.knots,...
                    nsub(rdim)-1, degree(rdim), regularity(rdim));
t_knots = kntunclamp(t_knots, degree(rdim), regularity(rdim), prdc_sides);

% Construct msh structures.
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

[xqn, xqw] = msh_set_quad_nodes (x_zeta, rule(1:end-1));
x_msh   = msh_cartesian (x_zeta, xqn, xqw, x_geo);

[tqn, tqw] = msh_set_quad_nodes (t_zeta, rule(end));
t_msh   = msh_cartesian (t_zeta, tqn, tqw, t_geo);

% Construct the space structures.
space   = sp_bspline (knots, degree, msh);
x_space = sp_bspline (x_knots, degree(1:rdim-1), x_msh);
t_space = sp_bspline (t_knots, degree(rdim), t_msh);

% Assembly the matrices.
Wt = op_gradu_v_tp (t_space, t_space, t_msh); %Controllo come è fatto!
Wt = Wt';
Mt = op_u_v_tp (t_space, t_space, t_msh);
Ms = op_u_v_tp (x_space, x_space, x_msh);
Ks = op_gradu_gradv_tp (x_space, x_space, x_msh);
A  = gmm*kron(Wt,Ms)+eta*kron(Mt,Ks);
R  = op_u_v_tp (space, space, msh);
S  = op_u_v_tp (space, space, msh);

% Extract internal dofs.
[~, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
int_dofs = setdiff (1:space.ndof, drchlt_dofs);
A = A(int_dofs,int_dofs);
R = R(int_dofs,int_dofs);
S = S(int_dofs,int_dofs);
L  = A'*(S\A);

% Compute the inf-sup 
[V, D] = eigs (L, R, 1, 'smallestabs');
mu = sqrt(real(D));
end