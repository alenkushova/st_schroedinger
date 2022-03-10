degmin = 3;
degmax = 7;
lmin = 3;
lmax = 7;

clear problem_data
geometry.time = geo_load(nrbline(0,1));
geometry.space = geo_load(nrbline(0,1));

for l = lmin:lmax
    
    fprintf('\n \n');
    fprintf('\n Discretization Level = %d \n ',l);
    fprintf('\n p     FD-backslash rel. err.     \n');
    fprintf('----------------------------------- \n');
     
    for deg = degmin:degmax

        degree_s     = deg;    % Degree of the splines
        regularity_s = deg-1;    % Regularity of the splines
        nsub_s       = 2^l;    % Number of subdivisions
        nquad_s      = deg+1;    % Points for the Gaussian quadrature rule
        degree_t     = deg;    % Degree of the splines
        regularity_t = deg-1;    % Regularity of the splines
        nsub_t       = 2^l;    % Number of subdivisions
        nquad_t      = deg+1;    % Points for the Gaussian quadrature rule
 
        [knots_s, zeta_s] = kntrefine (geometry.space.nurbs.knots, nsub_s-1, degree_s, regularity_s);
        [knots_t, zeta_t] = kntrefine (geometry.time.nurbs.knots, nsub_t-1, degree_t, regularity_t);
        
        % Construct msh structure
        [qn_s, qw_s] = msh_set_quad_nodes (zeta_s, msh_gauss_nodes (nquad_s));
        msh_s      = msh_cartesian (zeta_s, qn_s, qw_s, geometry.space);
        [qn_t, qw_t] = msh_set_quad_nodes (zeta_t, msh_gauss_nodes (nquad_t));
        msh_t      = msh_cartesian (zeta_t, qn_t, qw_t, geometry.time);
        
        % Construct space structure
        space_s    = sp_bspline (knots_s, degree_s, msh_s);
        space_t    = sp_bspline (knots_t, degree_t, msh_t);
        
        % Building of the SPACE/TIME matrices
        Ks = op_gradu_gradv_tp (space_s , space_s , msh_s, @(t)ones(size(t)));
        Ms = op_u_v_tp(space_s, space_s, msh_s);
        Wt = op_vel_dot_gradu_v_tp (space_t, space_t, msh_t, @(t)ones(size(t)));
        Mt = op_u_v_tp(space_t, space_t, msh_t);
        
        % simmetrizzo e impongo le b.c.
        Wt = Wt(2:end,2:end);
        Mt = (Mt(2:end,2:end) + Mt(2:end,2:end)')/2;
        Ms = (Ms(2:end-1,2:end-1) + Ms(2:end-1,2:end-1)')/2;
        Ks = (Ks(2:end-1,2:end-1) + Ks(2:end-1,2:end-1)')/2;
        
        nt = size(Wt,1);
        ns = size(Ks,1);
        F = randn(nt*ns,1);
              
        A = kron(Wt,Ms) + kron(Mt,Ks);
        u_bsl = A\F;

          % calcolo Ut con eig
%         [Ut,Lambda_t] = eig(full(Wt(1:end-1, 1:end-1)),full(Mt(1:end-1, 1:end-1)),'vector');
%         for j = 1:nt-1
%             Ut(:,j) = Ut(:,j)/sqrt(Ut(:,j)'*(Mt(1:end-1, 1:end-1)*Ut(:,j)));
%         end
        
        % calcolo Ut con schur
        L = chol(Mt(1:end-1, 1:end-1));
        [Ut, Lambda_t] = schur(full(((L')\Wt(1:end-1,1:end-1))/L),'complex');
        Lambda_t = sqrt(-1)*imag(diag(Lambda_t));
        Ut = L\Ut;
        
        v = -Mt(1:end-1,1:end-1)\Mt(1:end-1,end);
        r = [v;1]/sqrt([v' 1]*Mt*[v;1]);
        g = Ut'*Wt(1:end-1,:)*r;
        sigma = r'*Wt*r;
        Delta_t = sparse([diag(Lambda_t) g; -g' sigma]);
        Ut = [Ut; zeros(1,size(Ut,2))]; 
        Ut = [Ut r];
        
        % setup in space
        [Us,Lambda_s] = eig(full(Ks),full(Ms));
        for i =1:ns
            Us(:,i) = Us(:,i)/sqrt(Us(:,i)'*(Ms*Us(:,i))); 
        end
        Lambda_s = sparse(Lambda_s);
        
        H = kron(Delta_t,speye(ns)) + kron(speye(nt),Lambda_s);
        
        u_fd = Us'*reshape(F,ns,nt)*conj(Ut);
        u_fd = H\u_fd(:);
        u_fd = Us*reshape(u_fd,ns,nt)*Ut.';
        u_fd = u_fd(:);
        
        err = norm(u_fd - u_bsl)/norm(u_bsl);
        
        fprintf(' %d      %.3e \n', deg,err);
        
    end
end
        