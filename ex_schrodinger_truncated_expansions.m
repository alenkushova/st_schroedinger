%% EX_SCHRODINGER: solve the Schroedinger problem in one dimension. (S-T)
clear 
close all
clc
M = 3;

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% in this moment this value is fixed! Need to define better solver!
T = 2; %change T --> must modify 'geo_square_cilinder.txt' too !
%T = 10 %_______________________________________________________to use 'geo_square_cilinder.txt'
problem_data.T = T ;      % Final time.

% Physical domain, defined as NURBS map given in a text file
problem_data.xt_geo_name = 'geo_rectangle.txt'; %here final time is T = 2.
%problem_data.xt_geo_name = 'geo_square.txt'; %________to use 'geo_square_cilinder.txt'
problem_data.x_geo_name = nrbline ([0 0], [1 0]);
problem_data.t_geo_name = nrbline ([0 0], [T 0]);

% Type of boundary conditions for each side of the domain

problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 ];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 ];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Exact solution:
[uex, grad_uex, f] = Fourier_expansion(M);
%load('solutions'); % or you can load previously built solutions and rhs.

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
p = 3; %Polynomial degree
n = 256; % number of subdivisions in space direction!
method_data.degree     = [p p]; % Degree of the splines (last is time dir)
method_data.regularity = method_data.degree-1; % Regularity of the splines
method_data.nsub       = [n T*n]; % Number of subdivisions
method_data.nquad      = method_data.degree+1; % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';     % Fast Diag 'FD' or Matlab Backslash 'M'

%% 3) CALL TO THE SOLVER
tic
[geometry, msh, space, u] = solve_schrodinger_st_new (problem_data, method_data);
toc

%% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

output_file = 'Schroedinger_st_1D';

vtk_pts = {linspace(0, 1, 100), linspace(0, T, 100)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
subplot (1,2,1)
h1 = pcolor (X, Y, real(eu));
colorbar
colormap jet
h1.EdgeColor = 'none';
h1.FaceColor = 'interp';
title ('Numerical solution: \Re(u_h)'), axis tight
xlabel('x')
ylabel('Time')
subplot (1,2,2)
h2 = pcolor (X, Y, real(problem_data.uex (X,Y)));
colorbar
colormap jet
h2.EdgeColor = 'none';
h2.FaceColor = 'interp';
title ('Exact solution: \Re(u)'), axis tight
xlabel('x')
ylabel('Time')

%% 5) DISPLAY ERRORS of the computed solution in the L2 and H1 norm
% compute the error for the real part:
 Uex = @(x, t) real(problem_data.uex(x, t));
 GradUex = @(x, t) real(problem_data.graduex(x, t));
 [error_h1, error_l2] = sp_h1_error (space, msh, real(u), Uex, GradUex)

%% 6) Save solution
n = method_data.nsub;
d = method_data.degree;
filename = ['test_schrodinger_degree_' num2str(d) '_subs_' num2str(n) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);

