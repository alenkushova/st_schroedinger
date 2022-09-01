%% EX_SCHRODINGER: solve the Schroedinger problem in one dimension. Theta Scheme.
clear 
clc
M = 1250;

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% in this moment this value is fixed! Need to define better solver!
%T = 2; %change T --> must modify 'geo_square_cilinder.txt' too !
T = 1.e-3; % to be used with 'geo_square_cilinder.txt'


problem_data.T = T ;      % Final time.

% Physical domain, defined as NURBS map given in a text file
%problem_data.xt_geo_name = 'geo_rectangle.txt'; %here final time is T = 2.
problem_data.xt_geo_name = 'geo_small_rectangle.txt'; %here final time is T = 1.e-3.
%problem_data.xt_geo_name = 'geo_square.txt'; %________to use 'geo_square_cilinder.txt'
problem_data.x_geo_name = nrbline ([0 0], [1 0]);
problem_data.t_geo_name = nrbline ([0 0], [T 0]);

% Type of boundary conditions for each side of the domain

problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 ];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 ];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Exact solution:
% [uex, grad_uex, dudx, dudt, f] = Fourier_expansion(M);
% save('solutions1250.mat','uex', 'grad_uex', 'dudx', 'dudt', 'f')
solutions = ['solutions' num2str(M) '.mat'];
load(solutions); % or you can load previously built solutions and rhs.

problem_data.uex     = @(x, t) uex(x,t);
problem_data.graduex = @(x, t) grad_uex(x,t);

% Source term
problem_data.f = @(x, t) f(x,t);
% Dirichlet boundary conditions
problem_data.h = @(x, t, ind) problem_data.uex(x,t);
problem_data.x_h = @(x, ind) zeros (size (x)); % auxiliary function. 
problem_data.gmm = 1i;
problem_data.eta = 1;

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
p = 3;                                                                     % Polynomial degree
n = 8;                                                                   % number of subdivisions in space direction
N = 1024;         k = T/N;
method_data.N = N;                                                         % number of subdivisions in time direction
method_data.k = k;                                                         % mesh size in time direction
method_data.degree     = [p p];                                            % Degree of the splines (last is time dir)
method_data.regularity = method_data.degree-1;                             % Regularity of the splines
method_data.nsub       = [n N];                                            % Number of subdivisions
method_data.nquad      = method_data.degree+1;                             % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';                                             % Fast Diag 'FD' or Matlab Backslash 'M'
method_data.theta      = 1;                                              % Theta parameter for the theta scheme
%% 3) CALL TO THE SOLVER
fprintf ('Solving the problem with theta scheme... \n')
tic
[geometry, msh, space, theta_u] = solve_schrodinger_theta_scheme (problem_data, method_data);
toc

fprintf ('Solving the problem with DWF and FD... \n')
tic
[geo, mesh, spazio, u] = solve_schrodinger_dwf (problem_data, method_data);
toc

%compute the L2 projection of the solution
tic
fprintf ('Computing the mass matrix... \n')
M  = op_u_v_tp (spazio, spazio, mesh );
toc
fprintf ('Computing the rhs... \n')
rhs= op_f_v_tp (spazio, mesh, uex);
toc
fprintf ('Computing the L2 projection of the solution u... \n')
Pu = M\rhs;
toc

%% 5) DISPLAY ERRORS of the computed solution in the L2 and H1 norm
% compute the error for the real part:
fprintf ('Computing the errors of theta scheme in L2 norms... \n')
Uex = @(x, t) real(problem_data.uex(x, t));
GradUex = @(x, t) real(problem_data.graduex(x, t));
err_l2 = 0;
% errore in spazio
tic
for time = 0 : N
    % errore in norma L2 per soluzione p(x,t) per ogni istante t:
    error_l2 = sp_l2_error (space, msh, real(theta_u(:,time+1)), @(x) Uex(x,(time)*k));
    err_l2 = [err_l2 error_l2];
    fprintf('\n \n');
    fprintf('Errors in space at t_n with n = %d \n', (time)*k)
    toc
    fprintf('--------------------------------\n');
end
%errore in tempo: norma L infinito
errore_l2_INF = max(abs(err_l2));
%errore in tempo: norma L2
errore_l2_L2 = (k/2*sum((err_l2(1:end-1)).^2+(err_l2(2:end)).^2))^(1/2)    %this is what we need to check 
fprintf('\n \n');
toc

% compute the error for the real part:
fprintf ('Computing the error of solution in L2 norm... \n')
tic
u_error_l2 = sp_l2_error (spazio, mesh, real(u), Uex)
toc
% compute the error for the real part:
fprintf ('Computing the error of projection in L2 norm... \n')
tic
Pu_error_l2 = sp_l2_error (spazio, mesh, real(Pu), Uex)
toc

%% 6) Save solution
n = method_data.nsub(1);
d = method_data.degree(1);
theta = method_data.theta;
filename = ['test_schrodinger_degree_' num2str(d) '_subs_' num2str(n) '_and_' num2str(N) '_theta_' num2str(theta) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);

