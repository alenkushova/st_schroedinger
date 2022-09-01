% SOLVE_SCHRODINGER_ST_New: Solve a Schrodinger problem with a B-spline 
% discretization in a space-time isoparametric approach. This solver never 
% assembles the matrix associated with the problem, but exploits the tensor
% calculus of its Kronecker structure.
%
% The function solves the evolution problem
%
%     i d_t (u) - d_x( d_x (u)) = f    in Omega   = (0,1)^d x [0,T] 
%                             u = g    on Omega_0 = (0,1)^d x {0}
%                             u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_schrodinger_st_new (problem_data, method_data)
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
function [geometry, msh, space, u] = ...
              solve_schrodinger_dwf (problem_data, method_data) 
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
% in  order to avoid this matrix implementation:
% A = gmm*kron(Wt,Ms)+eta*kron(Mt,Ks);
% NN = size(Ws,1)^3*size(Wt,1);
% F = zeros(size(Ws,1)^3*size(Wt,1),1);
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

[~, x_drchlt_dofs] = sp_drchlt_l2_proj (x_space, x_msh, x_h, x_drchlt_sides);
x_int_dofs = setdiff (1:x_space.ndof, x_drchlt_dofs);

switch solver
    case 'FD'
        % With the preconditioner:
        Wt = Wt(2:end,2:end);% W - time
        Mt = Mt(2:end,2:end);% Mass in time
        Mt = (Mt+Mt')/2; %symmetrize
        Wt0= Wt(1:end-1,1:end-1);
        w  = Wt(1:end-1,end);
        %omega= Wt(end,end);
        Mt0= Mt(1:end-1,1:end-1);
        m  = Mt(1:end-1,end); 
        %mu = Mt(end,end);
        Ms = Ms(x_int_dofs,x_int_dofs); % Mass in space
        Ms = (Ms+Ms')/2; %symmetrize
        Ks = Ks(x_int_dofs,x_int_dofs); % Stiff in space  
        Ks = (Ks+Ks')/2; %symmetrize

        % Generalized diagonalizations:_____________________________ STEP 1
        [Us  , Ds]   = eig ( full(Ks)  , full(Ms)); 
        for i = 1 : size(Us,2)
            %normalize w.r.t. Ms
            Us(:,i) = Us(:,i)/sqrt(Us(:,i)' * Ms * Us(:,i)); 
        end 
        % decomposition in time using eig
%         [Ut0, Dt0] = eig ( full(Wt0), full(Mt0));
%         for i = 1 : size(Ut0,2)
%             %normalize w.r.t. Mt_0
%             Ut0(:,i) = Ut0(:,i)/sqrt(Ut0(:,i)' * Mt0 * Ut0(:,i)); 
%         end
        % decomposition in time using Schur
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

        %here is the arrow like decomposition
        Ut = [Ut0 r;0*r' rho];
        Lambda= sparse([diag(Dt0) g; -g' sigma]); 
        
        nt = size(Ut,1);
        Ns = size(Us,1);
        Ids = ones(Ns,1);
        ds  = diag(Ds);
        H   = gmm*diag(Lambda).' + eta*ds;
        B   = gmm*(g.').*Ids;
        S   = H(:,nt) - sum((((B').')./H(:,1:nt-1)).*B,2);
        Fnew = reshape(F(int_dofs),Ns,nt);
% apply the generalized eigenvector matrix _________________ STEP 1 
        tilde_Fnew = Us'*Fnew*((Ut').'); 
% solve the arrow like system ______________________________ STEP 2 
        tilde_u_nt = (-sum((((B').').*tilde_Fnew(:,1:nt-1)./H(:,1:nt-1)),2)...
                     + tilde_Fnew(:,nt))./S;
        tilde_u = [(tilde_Fnew(:,1:nt-1)-B.*tilde_u_nt)./H(:,1:nt-1) tilde_u_nt];
        %tilde_u = Arrow\tilde_Fnew(:); This works for sure
% apply the generalized eigenvector matrix _________________ STEP 3
        tilde_u = reshape(tilde_u,Ns,nt); %non dovrebbe servire.
%         u_sol = tmprod(tilde_u,{Ut Us},[2 1]); % may not be gpu Available
        u_sol = Us*tilde_u*(Ut.');
        u(int_dofs) = u_sol; % Numerical solution for Schrodinger.

    case 'M'
        A = gmm*kron(Wt,Ms)+eta*kron(Mt,Ks);
        u(int_dofs) = A(int_dofs,int_dofs)\F(int_dofs);
        
    otherwise
        A = gmm*kron(Wt,Ms)+eta*kron(Mt,Ks);
        u(int_dofs) = A(int_dofs,int_dofs)\F(int_dofs);
end

end