% this function sets the discrete tp B-spline spaces for the schroedinger
% problem, that is S^{p,p} space-time domain of maximum smoothness.
%
% Call: [msh, mshx, msht, space, x_space, t_space] = set_schrodinger_dwf (p,nel)
%
% Input: p   = polinomial degree
%        nel = number of elements per univariate direction
%
% Output:
%    geometry = geometry of the domain
%         msh = mesh associated to B-splines space-time 
%        mshx = mesh associated to B-splines in x direction(s)
%        msht = mesh associated to B-splines in time direction 
%       space = B-spline space 
%     x_space =   "        "   in x direction(s)
%     t_space =   "        "   in time direction
%

function [geometry, msh, mshx, msht, space, x_space, t_space] = set_schrodinger_dwf (p,nel)
% Geometries.
T = 1;
xt_geo_name = 'geo_square.txt';
x_geo_name  = nrbline ([0 0], [1 0]);
t_geo_name  = nrbline ([0 0], [T 0]);

geometry = geo_load(xt_geo_name);
geox    = geo_load(x_geo_name);
geot    = geo_load(t_geo_name);

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

[x_knots, x_zeta]= kntrefine(geox.nurbs.knots,...
                    nsub(1:rdim-1)-1, degree(1:rdim-1), regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, degree(1:rdim-1), regularity(1:rdim-1), prdc_sides);

[t_knots, t_zeta]= kntrefine(geot.nurbs.knots,...
                    nsub(rdim)-1, degree(rdim), regularity(rdim));
t_knots = kntunclamp(t_knots, degree(rdim), regularity(rdim), prdc_sides);

% Construct msh structures.
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

[xqn, xqw] = msh_set_quad_nodes (x_zeta, rule(1:end-1));
mshx   = msh_cartesian (x_zeta, xqn, xqw, geox);

[tqn, tqw] = msh_set_quad_nodes (t_zeta, rule(end));
msht   = msh_cartesian (t_zeta, tqn, tqw, geot);

% Construct the space structures.
space   = sp_bspline (knots, degree, msh);
x_space = sp_bspline (x_knots, degree(1:rdim-1), mshx);
t_space = sp_bspline (t_knots, degree(rdim), msht);

end