% This code is to compute the minimum eigenfunction w_min associated
% to the infsup condition of schroedinger operator with Graph norm in trial
% space and L2 norm in test space.
% 
%  p   =  degree of B-splines
%  nel =  number of subdivisions in space direction!
% 
% we check inf_v sup_w (A*v,w)/ |v||w| 
% A is the matrix associated to the bilinear form.
% Mv is the matrix representing the norm on trial space.
% Mw is the matrix representing the norm on test space.
% 
function [geo, msh, sp, res] = residual_rates(p, nel)
% Boundaries.
trial_drchlt_sides   = [1 2 3];% Dirichlet.
test_drchlt_sides   = [1 2 3];% Dirichlet.
h   = @(x, t, ind) 0*x;

% Set the discretization
input_args = set_schrodinger_dwf (nel, p, p);
dataset = set_discretization(input_args);

% Assembly the matrices.
A  = op_schroedinger_1st_order(dataset.xmsh,         dataset.tmsh,...
                               dataset.xspace_trial, dataset.tspace_trial,...
                               dataset.xspace_test,  dataset.tspace_test );

% trial Graph norm test L2
Mw = op_u_v_tp (dataset.xtspace_test, dataset.xtspace_test, dataset.xtmsh);
Mv = op_schroedinger_graph_norm (dataset.xtmsh, dataset.xmsh, dataset.tmsh,...
        dataset.xtspace_trial, dataset.xspace_trial, dataset.tspace_trial);

% symmetrize the mass matrices.
Mw = Mw+Mw'/2;    
Mv = Mv+Mv'/2;    

% Extract internal dofs.
[~, trial_drchlt_dofs] = sp_drchlt_l2_proj (dataset.xtspace_trial, dataset.xtmsh, h, trial_drchlt_sides);
trial_int_dofs = setdiff (1:dataset.xtspace_trial.ndof, trial_drchlt_dofs);
[~, test_drchlt_dofs] = sp_drchlt_l2_proj (dataset.xtspace_test, dataset.xtmsh, h, test_drchlt_sides);
test_int_dofs = setdiff (1:dataset.xtspace_test.ndof, test_drchlt_dofs);
Atilde  = A  (test_int_dofs, trial_int_dofs); %schroedinger op matrix
Mwtilde = Mw (test_int_dofs, test_int_dofs);  % norm in test space
Mvtilde = Mv (trial_int_dofs, trial_int_dofs);% norm in trial space

% Compute the inf-sup constants and eigenfunctions w
[mu, w, ~] = infsup_suriettivity(Atilde, Mvtilde, Mwtilde);
mu = sqrt(real(mu));

% compute v 
eigvect = zeros(dataset.xtspace_test.ndof,1);
v = zeros(dataset.xtspace_test.ndof,1);
eigvect(test_int_dofs,1) = w;
v(trial_int_dofs,1) = ((Mvtilde\(Atilde'))*w)/(mu^2);

% compute L2 projection on W of the field Sv
rhs = my_op_f_v_tp (dataset.xtspace_test, dataset.xtmsh, dataset.xtspace_trial, v);
Psv = Mw\rhs;

% compute ||(Id - P) Sv ||_L2
errl2 = schroedinger_l2_error_projection(dataset.xtspace_trial, dataset.xtmsh, v,...
                                 dataset.xtspace_test, Psv);
 
% compute || w ||_L2 
norm_w = abs(sqrt(eigvect.'*Mw'*eigvect));

% outputs 
res = errl2/norm_w;
geo = dataset.xtgeo;
msh = dataset.xtmsh;
sp  = dataset.xtspace_trial;

end
    
    
    
    
    
    
    
    
    
    
    
    
    
    