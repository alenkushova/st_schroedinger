%% EX_SCHRODINGER_B: solve the Schroedinger problem in one dimension. (S-T)
% THIS FILE HAVE PROBLEMS.
clear 
close all
clc

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% in this moment this value is fixed! Need to define better solver!
T = 1; %change T --> must modify 'geo_square_cilinder.txt' too !
%T = 10 %_______________________________________________________to use 'geo_square_cilinder.txt'
problem_data.T = T ;      % Final time.

% Physical domain, defined as NURBS map given in a text file
problem_data.xt_geo_name = 'geo_hypercube.txt'; %here final time is T = 1.
%problem_data.xt_geo_name = 'geo_square_cilinder.txt'; %________to use 'geo_square_cilinder.txt'
problem_data.x_geo_name = 'geo_cube.txt';
problem_data.t_geo_name = nrbline ([0 0], [T 0]);

% Type of boundary conditions for each side of the domain

problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 4 5 6 7];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 3 4 5 6];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Parameters:
omg = 0.2;
a0  = (2/(omg^2))^(1/4);

% Exact solution:
problem_data.uex     = @(x, y, z, t) (a0*exp(-1i*(x.^2 + y.^2 + z.^2 + t.^2)/(omg^2)));
problem_data.graduex = @(x, y, z, t) (cat (1, ...
                reshape (-2i*x/(omg^2).*problem_data.uex(x, y, z, t) , [1, size(x)]), ...
                reshape (-2i*y/(omg^2).*problem_data.uex(x, y, z, t) , [1, size(y)]), ...
                reshape (-2i*z/(omg^2).*problem_data.uex(x, y, z, t) , [1, size(z)]), ...
                reshape (-2i*t/(omg^2).*problem_data.uex(x, y, z, t) , [1, size(t)])));
% Source term
problem_data.f = @(x, y, z, t) ...
    zeros( size (x));
%            (4i*omg^2 + 4*x.^2 + 4*y.^2 + 2*omg^2*t)/omg^4.*problem_data.uex(x, y, t);
problem_data.h   = @(x, y, z, t, ind) problem_data.uex(x, y, z, t);
problem_data.x_h = @(x, y, z, ind) zeros (size (x)); % auxiliary function. 
problem_data.gmm = 1i;
problem_data.eta = 1;


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3 3 3]; % Degree of the splines (last is time dir)
method_data.regularity = method_data.degree-1; % Regularity of the splines
method_data.nsub       = [8 8 8 8]; % Number of subdivisions
method_data.nquad      = method_data.degree+1; % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';     % Fast Diag 'FD' or Matlab Backslash 'M'

%% 3) CALL TO THE SOLVER
fprintf('Solving the problem...')
tic
u = solve_schrodinger_dwf (problem_data, method_data);
toc
