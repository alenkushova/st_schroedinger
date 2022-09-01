function [x_geo, x_msh, x_space, u] = ...
              solve_schrodinger_theta_scheme (problem_data, method_data) 
% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end
n = numel(degree);

% Construct geometry structures
geometry = geo_load(xt_geo_name);
x_geo    = geo_load(x_geo_name);
t_geo    = geo_load(t_geo_name);

[knots, zeta]= kntrefine(geometry.nurbs.knots, nsub-1, degree, regularity);
knots = kntunclamp(knots, degree, regularity, prdc_sides);

[x_knots, x_zeta]= kntrefine(x_geo.nurbs.knots,...
                    nsub(1:n-1)-1, degree(1:n-1), regularity(1:n-1));
x_knots = kntunclamp(x_knots, degree(1:n-1), regularity(1:n-1), prdc_sides);

[t_knots, t_zeta]= kntrefine(t_geo.nurbs.knots,...
                    nsub(n)-1, degree(n), regularity(n));
t_knots = kntunclamp(t_knots, degree(n), regularity(n), prdc_sides);

% Construct msh structures
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

[xqn, xqw] = msh_set_quad_nodes (x_zeta, rule(1:end-1));
x_msh   = msh_cartesian (x_zeta, xqn, xqw, x_geo);

[tqn, tqw] = msh_set_quad_nodes (t_zeta, rule(end));
t_msh   = msh_cartesian (t_zeta, tqn, tqw, t_geo);

% Construct the space structures
space   = sp_bspline (knots, degree, msh);
x_space = sp_bspline (x_knots, degree(1:n-1), x_msh);
t_space = sp_bspline (t_knots, degree(n), t_msh);

% Apply Dirichlet bpoundary conditions + initial data.
u = zeros (x_space.ndof, N+1);
F = zeros (x_space.ndof, N+1); 
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides); % bisogna gestire meglio i dati al bordo...
u(drchlt_dofs(1:x_space.ndof),1) = u_drchlt(1:x_space.ndof);               % per ora facciamo che sono nulli (Dirichlet omogenei)
Ms = op_u_v_tp (x_space, x_space, x_msh);                                  % matrice di massa
Ks = op_gradu_gradv_tp (x_space, x_space, x_msh);                          % matrice di stiffness
for time = 0 : N 
    F(:,time+1) = op_f_v_tp (x_space, x_msh,@(x) f(x,time*k));             % here is the projection of f.
end
A = Ms + (k/1i)* (    theta) *Ks;
B = Ms - (k/1i)* (1 - theta) *Ks;

for time = 1 : N 
    rhs = B*u(:,time) + (k/1i)*(1-theta)*F(:,time) + (k/1i)*theta*F(:,time+1);
    u(2:end-1,time+1) = A(2:end-1,2:end-1)\rhs(2:end-1);
end

end
