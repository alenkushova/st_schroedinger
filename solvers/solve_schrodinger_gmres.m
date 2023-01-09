% SOLVE_SCHRODINGER_GMRES: Solve a Schrodinger problem with a B-spline 
% discretization in a space-time isoparametric approach. This solver  
% assembles the matrix associated with the problem, and uses an iterative
% gmres solver preconditioned with the kroneker matrix decomposed as arrow
% like extended fast diagonalization. 
%
% The function solves the evolution problem
%
%     i d_t (u) - d_x( d_x (u)) = f    in Q   = Omega^d x [0,T] 
%                             u = g    on Q_0 = Omega^d x {0}
%                             u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_schrodinger_gmres (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition 
%    - prdc_sides:   sides with Periodic boundary condition (may be empty)
%    - f:            function for the source term (right hand side)
%    - g:            function for Neumann boundary condition (if needed)
%    - h:            function for Dirichlet boundary condition (also initial)
%    - T:            final time T.
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
function [geometry, msh, space, u, gmres_output] = ...
              solve_schrodinger_gmres (problem_data, method_data) 
          
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

% Assembly the matrices
Wt = op_gradu_v_tp (t_space, t_space, t_msh); %Controllo come è fatto!
Wt = Wt';
Mt = op_u_v_tp (t_space, t_space, t_msh);
Ms = op_u_v_tp (x_space, x_space, x_msh);
Ks = op_gradu_gradv_tp (x_space, x_space, x_msh);

%A = gmm*kron(Wt,Ms)+eta*kron(Mt,Ks); %the matrix associated to the problem.
% with 128^3 elements are Requested 2248091x2248091 (16.6GB) 
% the array exceeds maximum array size preference (15.7GB). 

fprintf('building rhs... \n\n')
%right hand side
F = op_f_v_tp (space, msh, f); %here is the projection of f.

% Apply Dirichlet bpoundary conditions
u = zeros (space.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;
int_dofs = setdiff (1:space.ndof, drchlt_dofs);

% we modify the rhs this way:
mat_u = reshape(u,x_space.ndof,t_space.ndof);
v1 = gmm*Ms*mat_u*(Wt.') + eta*Ks*mat_u*(Mt.');
v1 = v1(:);
F(int_dofs) = F(int_dofs) - v1(int_dofs); 
% Solution with backslash
% u(int_dofs) = A(int_dofs,int_dofs)\F(int_dofs);
% ubs = u;

[~, x_drchlt_dofs] = sp_drchlt_l2_proj (x_space, x_msh, x_h, x_drchlt_sides);
x_int_dofs = setdiff (1:x_space.ndof, x_drchlt_dofs);

% % Build explicitly the preconditioner:
% input_args = set_input_arguments (n, degree(1), degree(1), nsub(1));
% dataset    = set_discretization(input_args);
% P  = set_preconditioner (dataset.xmsh, dataset.tmsh, ...
%       dataset.xspace_trial, dataset.tspace_trial);

% Build decomposition U(arrow)U* of the preconditioner:
input_args = set_input_arguments (2, degree(1), degree(1), nsub(1));
dataset    = set_univ_spaces_xt(input_args);
[Ut, Us, dt, g, ds, dtg] = GenFastDiag(dataset.xmsh, dataset.tmsh, ...
                                        dataset.xspace, dataset.tspace) ;
%[Ut, Us, dt, g, ds, dtg] = GenFastDiag (x_msh, t_msh, x_space, t_space, x_space, t_space);
                                    
nt = numel(dt);
ns = numel(ds)^(numel(degree)-1);
ms = Ms(x_int_dofs,x_int_dofs);
mt = Mt(2:end,2:end);
wt = Wt(2:end,2:end);
ks = Ks(x_int_dofs,x_int_dofs);
b = F(int_dofs);
%c = @(x) A(int_dofs,int_dofs)\x; %con questa ovviamente ci metto una
%iterazione soltanto. Dovrei scrivere bene la funzione sotto SolveFD
tol = 10^(-8);
maxit = min(size(b,1),1000);
switch n
  case 2
    a = @(x) op_Ax1D(ms,mt,wt,ks,ns,nt,gmm,eta,x); %questa funzione applica A ed è corretta
    c =@(x) SolveFD_1D(Us,Ut,dt,g,ds,gmm,eta,x);
  case 3
    a = @(x) op_Ax2D(ms,mt,wt,ks,ns,nt,gmm,eta,x); %questa funzione applica A ed è corretta
    c =@(x) SolveFD_2D(Us,Ut,dt,g,ds,gmm,eta,x);
  case 4 
    a = @(x) op_Ax3D(ms,mt,wt,ks,ns,nt,gmm,eta,x); %questa funzione applica A ed è corretta
    c =@(x) SolveFD_3D(Us,Ut,dt,g,ds,gmm,eta,x);
end
% Ds = kron(diag(ds),speye(sqrt(ns))) + kron(speye(sqrt(ns)),diag(ds));
% arrow = gmm*kron(dtg,speye(ns)) + eta*kron(speye(nt),Ds);
% Utrasp = kron(Ut',kron(Us',Us'));
% U = kron(Ut,kron(Us,Us));
% c = @(x) Utrasp*(arrow\(U*x));
fprintf('solving... \n\n')
[u_inner, flag, rel_res, iter, res_vec] = gmres(a,b,[],tol,maxit,c);
u(int_dofs) = u_inner;
% Output informations of GMRES solver iterations
gmres_output.flag    = flag;
gmres_output.rel_res = rel_res;
gmres_output.iter    = iter;
gmres_output.res_vec = res_vec;

end

function input_args = set_input_arguments (n, trial_p, test_p, nsub)
switch n
    case 2
        input_args.xt_geo_name = 'geo_square.txt';
        input_args.x_geo_name  = nrbline ([0 0], [1 0]);
    case 3
        input_args.xt_geo_name = 'geo_cube.txt';        
        input_args.x_geo_name  = 'geo_square.txt';
    case 4
        input_args.xt_geo_name = 'geo_hypercube.txt';
        input_args.x_geo_name  = 'geo_cube.txt';
end
input_args.t_geo_name  = nrbline ([0 0], [1 0]); % Only T = 1 implemented yet.
input_args.trial_degree     = repmat(trial_p, 1, n);                       % Degree of the splines (last is time dir)
input_args.trial_regularity = input_args.trial_degree-1;                   % Regularity of the splines ( " )
input_args.test_degree      = repmat(test_p, 1, n);                        % Degree of the splines (last is time dir)
input_args.test_regularity  = input_args.test_degree-1;                    % Regularity of the splines ( " )
input_args.nsub  = repmat(nsub, 1, n);                                     % Number of subdivisions ( " )
input_args.nquad = max(input_args.trial_degree, input_args.test_degree) +1;% Points for the Gaussian quadrature rule ( " )
end

function dataset    = set_discretization(input_args)

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

% space time knots for test functions
[knots, ~]    = kntrefine(dataset.xtgeo.nurbs.knots, nsub-1, test_degree, test_regularity);
knots   = kntunclamp(knots, test_degree, test_regularity, []);

% space knots for test functions in space
[x_knots, ~]= kntrefine(dataset.xgeo.nurbs.knots,...
                    nsub(1:rdim-1)-1, test_degree(1:rdim-1), test_regularity(1:rdim-1));
x_knots = kntunclamp(x_knots, test_degree(1:rdim-1), test_regularity(1:rdim-1), []);

%time knots for test functions in time
[t_knots, ~]= kntrefine(dataset.tgeo.nurbs.knots,...
                    nsub(rdim)-1, test_degree(rdim), test_regularity(rdim));
t_knots = kntunclamp(t_knots, test_degree(rdim), test_regularity(rdim), []);

% define space structures for test functions
dataset.xtspace_test = sp_bspline (knots, test_degree, dataset.xtmsh);
dataset.xspace_test  = sp_bspline (x_knots, test_degree(1:rdim-1), dataset.xmsh);
dataset.tspace_test  = sp_bspline (t_knots, test_degree(rdim), dataset.tmsh);

end

function varargout  = set_preconditioner (mshx, msht, trial_spx, trial_spt,...
                                           test_spx, test_spt)
      
  if (nargin < 4)
    error('op_schroedinger_1st_order: not enought imput arguments');
  elseif (nargin == 4)
    test_spx = trial_spx;
    test_spt = trial_spt;
  elseif (nargin == 5 ) 
    test_spt = trial_spt;
  end
      
  Wt = op_gradu_v_tp (trial_spt, test_spt, msht);
  Mt = op_u_v_tp (trial_spt, test_spt, msht);
  Ms = op_u_v_tp (trial_spx, test_spx, mshx);
  Ks = op_gradu_gradv_tp (trial_spx, test_spx, mshx);
  A  = 1i*kron(Wt',Ms) + kron(Mt,Ks);
  
  if (nargout == 1 || nargout == 0)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, values] = find(A);
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_schroedinger_1st_order: wrong number of output arguments')
  end

end

function dataset = set_univ_spaces_xt(input_args)
% read the fields from structures
data_names = fieldnames (input_args);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= input_args.(data_names{iopt});']);
end

dataset.xgeo  = geo_load(x_geo_name);
dataset.tgeo  = geo_load(t_geo_name);

% define mesh structures :
rdim = numel(trial_degree);

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

[xqn, xqw]   = msh_set_quad_nodes (x_zeta, rule(1:end-1));
dataset.xmsh = msh_cartesian (x_zeta, xqn, xqw, dataset.xgeo);

[tqn, tqw]   = msh_set_quad_nodes (t_zeta, rule(end));
dataset.tmsh = msh_cartesian (t_zeta, tqn, tqw, dataset.tgeo);

% define space structures for trial functions
dataset.xspace  = sp_bspline (x_knots, trial_degree(1:rdim-1), dataset.xmsh);
dataset.tspace  = sp_bspline (t_knots, trial_degree(rdim), dataset.tmsh);

end

function [Ut, Us, dt, g, ds, dtg] = GenFastDiag (mshx, msht, ... 
                                  trial_spx, trial_spt, test_spx, test_spt)
      
  if (nargin < 4)
    error('op_schroedinger_1st_order: not enought imput arguments');
  elseif (nargin == 4)
    test_spx = trial_spx;
    test_spt = trial_spt;
  elseif (nargin == 5 ) 
    test_spt = trial_spt;
  end
  
  Wt = op_gradu_v_tp (trial_spt, test_spt, msht); Wt = Wt';
  Mt = op_u_v_tp (trial_spt, test_spt, msht);
  Ms = op_u_v_tp (trial_spx, test_spx, mshx);
  Ks = op_gradu_gradv_tp (trial_spx, test_spx, mshx);
  
  % GENERALIZED FAST DIAGONALIZATION: 
  % IN TIME...
  Wt = Wt(2:end,2:end);    % W in time
  Mt = Mt(2:end,2:end);    % Mass in time
  Mt = (Mt+Mt')/2;         % Symmetrize
  Wt0= Wt(1:end-1,1:end-1);
  w  = Wt(1:end-1,end);
  %omega= Wt(end,end);
  Mt0= Mt(1:end-1,1:end-1);
  m  = Mt(1:end-1,end);
  %mu = Mt(end,end);
  %IN SPACE...
  Ms = Ms(2:end-1,2:end-1); % Mass in space
  Ms = (Ms+Ms')/2;          % Symmetrize
  Ks = Ks(2:end-1,2:end-1); % Stiff in space
  Ks = (Ks+Ks')/2;          % Symmetrize
  % generalized eigendecomposition in space STEP 1
  [Us , Ds]   = eig ( full(Ks)  , full(Ms));
  ds = diag(sparse(Ds));
  for i = 1 : size(Us,2)
      Us(:,i) = Us(:,i)/sqrt(Us(:,i)' * Ms * Us(:,i)); %normalize w.r.t. Ms
  end
  % generalized eigendecomposition in time using Schur STEP 2
  L = chol(Mt0);
  [Ut0, Dt0] = schur(full(((L')\Wt0)/L),'complex');
  Dt0 = sqrt(-1)*imag(diag(Dt0));
  Ut0 = L\Ut0;
  
  %find v
  v = Mt0 \ (-m);
  v_1   = [v;1];

  %find r and rho
  r_rho = v_1/sqrt(v_1' * Mt * v_1);
  r     = r_rho(1:end-1);
  rho   = r_rho(end);
  
  %find g and sigma
  g     = Ut0' * [Wt0 w] *r_rho;
  sigma = r_rho' * Wt * r_rho;

  % here is the arrow like decomposition
  Ut = [Ut0 r;0*r' rho];
  dtg= sparse([diag(Dt0) g; -g' sigma]);
  dt = diag(dtg);
end

function u = SolveFD_1D(Us,Ut,dt,g,ds,gmm,eta,rhs)

ns = size(Us,1); nt = size(Ut,1); Ns = ns;                               % taking the sizes
                  
Ds = ds;                                            % vector kronecker product of space eigenvalues

H  = gmm*dt.'+eta*Ds;                                                      % assembling the block matrices H
B  = gmm*(g.').*ones(Ns,1);                                                % assembling the block matrices B
S  = H(:,nt) - sum((((B').')./H(:,1:nt-1)).*B,2);                          % assembling the block matrix S

rhs = reshape(rhs, ns, nt);

%STEP 1
tilde_rhs = tmprod(rhs, {Us' Ut'}, 1:2);  
%tilde_rhs = Us'*rhs*((Ut').');
tilde_rhs = reshape (tilde_rhs, Ns, nt);
%STEP 2 - potrebbe necessitare modifiche...
tilde_u     = (-sum((((B').').*tilde_rhs(:,1:nt-1)./H(:,1:nt-1)),2)...
               + tilde_rhs(:,nt))./S;
new_tilde_u = [(tilde_rhs(:,1:nt-1)-B.*tilde_u)./H(:,1:nt-1) tilde_u];

%STEP 3
new_tilde_u = reshape(new_tilde_u, ns, nt);
%u = Us*new_tilde_u*(Ut.');
u = tmprod(new_tilde_u, {Us Ut}, 1:2);  
u=u(:);

end

function u = SolveFD_2D(Us,Ut,dt,g,ds,gmm,eta,rhs)

ns = size(Us,1); nt = size(Ut,1); Ns = ns^2;                               % taking the sizes
                  
if nargin < 8
    rhs = ones(ns,ns,nt);                                                  % if the right hand side is not given it is set to be ones
else
    rhs = reshape(rhs, ns, ns, nt);
end

id = ones(size(ds));
Ds = kron(ds,id) + kron(id,ds);                                            % vector kronecker product of space eigenvalues

H  = gmm*dt.'+eta*Ds;                                                      % assembling the block matrices H
B  = gmm*(g.').*ones(Ns,1);                                                % assembling the block matrices B
S  = H(:,nt) - sum((((B').')./H(:,1:nt-1)).*B,2);                          % assembling the block matrix S

%STEP 1
tilde_rhs = tmprod(rhs, {Us' Us' Ut'}, 1:3);  
tilde_rhs = reshape (tilde_rhs, Ns, nt);

%STEP 2 - potrebbe necessitare modifiche...
tilde_u     = (-sum((((B').').*tilde_rhs(:,1:nt-1)./H(:,1:nt-1)),2)...
               + tilde_rhs(:,nt))./S;
new_tilde_u = [(tilde_rhs(:,1:nt-1)-B.*tilde_u)./H(:,1:nt-1) tilde_u];

%STEP 3
new_tilde_u = reshape(new_tilde_u, ns, ns, nt);
u = tmprod(new_tilde_u, {Us Us Ut}, 1:3);  
u = u(:);

end

function u = SolveFD_3D(Us,Ut,dt,g,ds,gmm,eta,rhs)

ns = size(Us,1); nt = size(Ut,1); Ns = ns^3;                               % taking the sizes
                  
if nargin < 8
    rhs = ones(ns, ns, ns, nt);                                                  % if the right hand side is not given it is set to be ones
else
    rhs = reshape(rhs, ns, ns, ns, nt);
end

id = ones(size(ds));
Ds = kron(kron(ds,id),id) + kron(kron(id,ds),id) + kron(kron(id,id),ds);   % vector kronecker product of space eigenvalues

H  = gmm*dt.'+eta*Ds;                                                      % assembling the block matrices H
B  = gmm*(g.').*ones(Ns,1);                                                % assembling the block matrices B
S  = H(:,nt) - sum((((B').')./H(:,1:nt-1)).*B,2);                          % assembling the block matrix S

%STEP 1
tilde_rhs = tmprod(rhs, {Us' Us' Us' Ut'}, 1:4);  
tilde_rhs = reshape (tilde_rhs, Ns, nt);

%STEP 2 - potrebbe necessitare modifiche...
tilde_u     = (-sum((((B').').*tilde_rhs(:,1:nt-1)./H(:,1:nt-1)),2)...
               + tilde_rhs(:,nt))./S;
new_tilde_u = [(tilde_rhs(:,1:nt-1)-B.*tilde_u)./H(:,1:nt-1) tilde_u];

%STEP 3
new_tilde_u = reshape(new_tilde_u, ns, ns, ns, nt);
u = tmprod(new_tilde_u, {Us Us Us Ut}, 1:4);  
u = u(:);

end

function y = op_Ax1D(ms,mt,wt,ks,ns,nt,gmm,eta,x)
x = reshape(x,ns,nt);
y = gmm*ms*x*(wt.')+eta*ks*x*(mt.');
y = y(:);
end

function y = op_Ax2D(ms,mt,wt,ks,ns,nt,gmm,eta,x)
x = reshape(x,ns,nt);
y = gmm*ms*x*(wt.')+eta*ks*x*(mt.');
y = y(:);
end

function y = op_Ax3D(ms,mt,wt,ks,ns,nt,gmm,eta,x)
x = reshape(x,ns,nt);
y = gmm*ms*x*(wt.')+eta*ks*x*(mt.');
y = y(:);
end

