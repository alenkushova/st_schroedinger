% This function sets the discrete tp B-spline spaces.
% 
% Call: dataset = set_discretization (input_args)
%
% Input: input_args is a structure with the following fields
%
% Output: 'dataset' structure with the following fields
%
%      geometry = geometry of the domain
%           msh = mesh associated to B-splines space-time 
%          mshx = mesh associated to B-splines in x direction(s)
%          msht = mesh associated to B-splines in time direction 
% xtspace_trial = B-spline space for trial functions
% x_space_trial =   "        "   in x direction(s)
% t_space_trial =   "        "   in time direction
%  xtspace_test = B-spline space for test functions
%  x_space_test =   "        "   in x direction(s)
%  t_space_test =   "        "   in time direction
%

function dataset = set_discretization(input_args)

% read the fields from structures
data_names = fieldnames (input_args);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= input_args.(data_names{iopt});']);
end

% define geometies
dataset.xtgeo = geo_load(xt_geo_name);
dataset.xgeo  = geo_load(x_geo_name);
dataset.tgeo  = geo_load(t_geo_name);

% define mesh structures :
rdim = numel(trial_degree);
% space time knots for trial functions
[knots, zeta]    = kntrefine(dataset.xtgeo.nurbs.knots, nsub-1, trial_degree, trial_regularity);
knots   = kntunclamp(knots, trial_degree, trial_regularity, []);

% space knots for trial functions in space
[x_knots, x_zeta]= kntrefine(dataset.xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, trial_degree(1:rdim-1), trial_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, trial_degree(1:rdim-1), trial_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, t_zeta]= kntrefine(dataset.tgeo.nurbs.knots,...
                    nsub(rdim)-1, trial_degree(rdim), trial_regularity(rdim));
t_knots = kntunclamp(t_knots, trial_degree(rdim), trial_regularity(rdim), []);

% define quadrature rules 
rule         = msh_gauss_nodes (nquad);
[qn, qw]     = msh_set_quad_nodes (zeta, rule);
dataset.xtmsh= msh_cartesian (zeta, qn, qw, dataset.xtgeo);

[xqn, xqw]   = msh_set_quad_nodes (x_zeta, rule(1:end-1));
dataset.xmsh = msh_cartesian (x_zeta, xqn, xqw, dataset.xgeo);

[tqn, tqw]   = msh_set_quad_nodes (t_zeta, rule(end));
dataset.tmsh = msh_cartesian (t_zeta, tqn, tqw, dataset.tgeo);

% define space structures for trial functions
dataset.xtspace_trial = sp_bspline (knots, trial_degree, dataset.xtmsh);
dataset.xspace_trial  = sp_bspline (x_knots, trial_degree(1:rdim-1), dataset.xmsh);
dataset.tspace_trial  = sp_bspline (t_knots, trial_degree(rdim), dataset.tmsh);

% space time knots for trial functions
[knots, ~]    = kntrefine(dataset.xtgeo.nurbs.knots, nsub-1, test_degree, test_regularity);
knots   = kntunclamp(knots, test_degree, test_regularity, []);

% space knots for trial functions in space
[x_knots, ~]= kntrefine(dataset.xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, test_degree(1:rdim-1), test_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, test_degree(1:rdim-1), test_regularity(1:rdim-1), []);

%time knots for trial functions in time
[t_knots, ~]= kntrefine(dataset.tgeo.nurbs.knots,...
                    nsub(rdim)-1, test_degree(rdim), test_regularity(rdim));
t_knots = kntunclamp(t_knots, test_degree(rdim), test_regularity(rdim), []);

% define space structures for trial functions
dataset.xtspace_test = sp_bspline (knots, test_degree, dataset.xtmsh);
dataset.xspace_test  = sp_bspline (x_knots, test_degree(1:rdim-1), dataset.xmsh);
dataset.tspace_test  = sp_bspline (t_knots, test_degree(rdim), dataset.tmsh);

end