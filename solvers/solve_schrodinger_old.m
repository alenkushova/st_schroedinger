% SOLVE_SCHRODINGER_ST: Solve a Schrodinger problem with a B-spline
%                               discretization in a space-time isoparametric
%                               approach . 
%
% The function solves the evolution problem
%
%     i d_t (u) - d_x( d_x (u)) = f    in Omega = F((0,1)^(n+1)) 
%                             u = g    on Omega_0 = F((0,1)^n x {0})
%                             u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_schrodinger_st (problem_data, method_data)
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
%    - dimension:    dimension of the space variable (JUST 1 FOR NOW!!!!!!!!!!!!!!!!!!!!!!!!!!)
%    - T:            final time
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
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
function [GEO, MSH, SPACE, u] = ...
              solve_schrodinger_old (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

%Construct 1d geometry structures
% Construct geometry
geometry = geo_load(geo_name);

[knots, zeta]= kntrefine(geometry.nurbs.knots, nsub-1, degree, regularity);
knots = kntunclamp(knots, degree, regularity, prdc_sides);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
msh      = msh_cartesian (zeta, qn, qw, geometry);

% Construct space structure
space_t = sp_bspline (knots, degree, msh);
space_s = sp_bspline (knots, degree, msh);

% Tensorize all one dimentional spaces: (modify to deal with higher dim)
% NMNN_sides   = []; % Global Neumann 
DRCHLT_sides = [1 2 3];  % Global Dirichlet

% Construct geometry
DEGREE = degree      * ones(1,dimension + 1);
REG    = regularity  * ones(1,dimension + 1);
NSUB   = nsub        * ones(1,dimension + 1);
NQUAD  = nquad       * ones(1,dimension + 1);

GEO    = geo_load('geo_square.txt');

[KNOTS, ZETA]= kntrefine(GEO.nurbs.knots, NSUB-1, DEGREE, REG);
KNOTS  = kntunclamp(KNOTS, DEGREE, REG, []); 

% Construct global msh structure
RULE     = msh_gauss_nodes (NQUAD);
[QN, QW] = msh_set_quad_nodes (ZETA, RULE);
MSH      = msh_cartesian (ZETA, QN, QW, GEO);

% Construct global space structure
SPACE = sp_bspline (KNOTS, DEGREE, MSH);

% Assembly the matrices
Wt = op_gradu_v_tp     (space_t, space_t, msh); %Controllo come è fatto!
Wt = Wt';
Mt = op_u_v_tp         (space_t,space_t,msh);
Mt = T*Mt; 
Ms = op_u_v_tp         (space_s,space_s,msh);
Ks = op_gradu_gradv_tp (space_s, space_s, msh);
A = 1i*kron(Wt,Ms)+kron(Mt,Ks);
F = op_f_v_tp (SPACE, MSH, f); %here is the projection of f.

% Apply Dirichlet bpoundary conditions
u = zeros (SPACE.ndof, 1);
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (SPACE, MSH, h, DRCHLT_sides);
u(drchlt_dofs) = u_drchlt;
int_dofs = setdiff (1:SPACE.ndof, drchlt_dofs);
F(int_dofs) = F(int_dofs) - A(int_dofs, drchlt_dofs)*u_drchlt; %modify rhs.

switch solver
    case 'FD'
        % With the preconditioner:
        Wt   = Wt(2:end,2:end);% W - time
        Mt   = Mt(2:end,2:end);% Mass in time
        Mt = (Mt+Mt')/2; %symmetrize
        Wt0  = Wt(1:end-1,1:end-1);
        w    = Wt(1:end-1,end);
        %omega= Wt(end,end);
        Mt0  = Mt(1:end-1,1:end-1);
        m    = Mt(1:end-1,end); 
        %mu   = Mt(end,end);
        Ms   = Ms(2:end-1,2:end-1); % Stiff in space
        Ms = (Ms+Ms')/2; %symmetrize
        Ks = Ks(2:end-1,2:end-1); % Mass in space  
        Ks = (Ks+Ks')/2; %symmetrize

        % Generalized diagonalizations:
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
        U     = kron(Ut,Us);
        Arrow = 1i*kron(Lambda,eye(size(Ds,1)))+kron(eye(size(Lambda,1)),Ds);

        % Solve the system:
        tilde_F = U'*F(int_dofs); %STEP 2
        tilde_u = Arrow\tilde_F; % STEP 3 
        u(int_dofs) = U*tilde_u; % STEP 4 - Num. solution for Schrodinger
    case 'M'
        u(int_dofs) = A(int_dofs,int_dofs)\F(int_dofs);
    otherwise
        u(int_dofs) = A(int_dofs,int_dofs)\F(int_dofs);
end

end