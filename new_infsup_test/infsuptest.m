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
function [mu, D, geo, msh, space, eigvect,lambda2_Mw_w, v, Av] = ...
                                                       infsuptest (p, nel)
% Boundaries.
trial_drchlt_sides   = [1 2 3];% Dirichlet.
test_drchlt_sides   = [1 2 3];% Dirichlet.
h   = @(x, t, ind) 0*x;

% Set the discretization
input_args = set_schrodinger_dwf (nel, p, p);
%input_args = set_schrodinger_dsf (nel, p, p);
dataset = set_discretization(input_args);

% Assembly the matrices.
A  = op_schroedinger_1st_order(dataset.xmsh,         dataset.tmsh,...
                               dataset.xspace_trial, dataset.tspace_trial,...
                               dataset.xspace_test,  dataset.tspace_test );
% trial L2 norm test L2
Mw = op_u_v_tp (dataset.xtspace_test, dataset.xtspace_test, dataset.xtmsh);
Mv = op_u_v_tp (dataset.xtspace_trial, dataset.xtspace_trial, dataset.xtmsh);

% trial Graph norm test L2
% Mw = op_u_v_tp (dataset.xtspace_test, dataset.xtspace_test, dataset.xtmsh);
% Mv = op_schroedinger_graph_norm (dataset.xtmsh, dataset.xmsh, dataset.tmsh,...
%         dataset.xtspace_trial, dataset.xspace_trial, dataset.tspace_trial);

% trial L2 and test Graph norm
% Mv = op_u_v_tp (dataset.xtspace_trial, dataset.xtspace_trial, dataset.xtmsh);
% Mw = op_schroedinger_graph_norm (dataset.xtmsh, dataset.xmsh, dataset.tmsh,...
%        dataset.xtspace_test, dataset.xspace_test, dataset.tspace_test);

% Extract internal dofs.
[~, trial_drchlt_dofs] = sp_drchlt_l2_proj (dataset.xtspace_trial, dataset.xtmsh, h, trial_drchlt_sides);
trial_int_dofs = setdiff (1:dataset.xtspace_trial.ndof, trial_drchlt_dofs);
[~, test_drchlt_dofs] = sp_drchlt_l2_proj (dataset.xtspace_test, dataset.xtmsh, h, test_drchlt_sides);
test_int_dofs = setdiff (1:dataset.xtspace_test.ndof, test_drchlt_dofs);
Atilde  = A  (test_int_dofs, trial_int_dofs); %schroedinger op matrix
Mwtilde = Mw (test_int_dofs, test_int_dofs);  % norm in test space
Mvtilde = Mv (trial_int_dofs, trial_int_dofs);% norm in trial space

% Compute the inf-sup constants and eigenfunctions v:
% [mu, vtilde, D] = infsup_coercivity(Atilde, Mvtilde, Mwtilde);
% mu = sqrt(real(mu));
% eigvect = zeros(dataset.xtspace_test.ndof,1);
% v = zeros(dataset.xtspace_trial.ndof,1);
% Av = zeros(dataset.xtspace_test.ndof,1);
% v(trial_int_dofs,1) = vtilde;
% eigvect(test_int_dofs,1) = Mwtilde\ (Atilde*vtilde);

% Compute the inf-sup constants and eigenfunctions w
[mu, w, D] = infsup_suriettivity(Atilde, Mvtilde, Mwtilde);
mu = sqrt(real(mu));
eigvect = zeros(dataset.xtspace_test.ndof,1);
v = zeros(dataset.xtspace_test.ndof,1);
Av = zeros(dataset.xtspace_test.ndof,1);
eigvect(test_int_dofs,1) = w;
v(trial_int_dofs,1) = (Mvtilde\(Atilde'))*w;
Av(test_int_dofs,1) = (Atilde*v(trial_int_dofs,1));
lambda2_Mw_w = zeros(dataset.xtspace_test.ndof,1); 
lambda2_Mw_w(test_int_dofs,1) = mu*Mwtilde*w;

geo = dataset.xtgeo;
msh = dataset.xtmsh;
space = dataset.xtspace_trial;

end