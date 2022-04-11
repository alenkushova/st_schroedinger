%% EX_SCHRODINGER_B: solve the Schroedinger problem in one dimension. (S-T)
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
problem_data.xt_geo_name = 'geo_cube.txt'; %here final time is T = 1.
%problem_data.xt_geo_name = 'geo_square_cilinder.txt'; %________to use 'geo_square_cilinder.txt'
problem_data.x_geo_name = 'geo_square.txt';
problem_data.t_geo_name = nrbline ([0 0], [T 0]);

% Type of boundary conditions for each side of the domain

problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 4 5];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 3 4];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Parameters:
omg = 0.2;
a0  = (2/(omg^2))^(1/4);

% Exact solution:
problem_data.uex     = @(x, y, t) (a0*exp(-1i*(x.^2 + y.^2 + t.^2)/(omg^2)));
problem_data.graduex = @(x, y, t) (cat (1, ...
                reshape (-2i*x/(omg^2).*problem_data.uex(x, y, t) , [1, size(x)]), ...
                reshape (-2i*y/(omg^2).*problem_data.uex(x, y, t) , [1, size(y)]), ...
                reshape (-2i*t/(omg^2).*problem_data.uex(x, y, t) , [1, size(t)])));
% Source term
problem_data.f = @(x, y, t) ...
            (4i*omg^2 + 4*x.^2 + 4*y.^2 + 2*omg^2*t)/omg^4.*problem_data.uex(x, y, t);
problem_data.h = @(x, y, t, ind) problem_data.uex(x, y, t);
problem_data.x_h = @(x, y, ind) zeros (size (x)); % auxiliary function. 
problem_data.gmm = 1i;
problem_data.eta = 1;


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3 3]; % Degree of the splines (last is time dir)
method_data.regularity = method_data.degree-1; % Regularity of the splines
method_data.nsub       = [8 8 8]; % Number of subdivisions
method_data.nquad      = method_data.degree+1 ; % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';     % Fast Diag 'FD' or Matlab Backslash 'M'

%% 3) CALL TO THE SOLVER
tic
[geometry, msh, space, u] = solve_schrodinger_st_new (problem_data, method_data);
toc
%% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW

% output_file = 'Schroedinger_st_2D';

vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

figure
sp_plot_solution (real(u), space, geometry, [40 40 40], [20 20 20]);

% %% 4) POSTPROCESSING
% vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)};
% [eu, F] = sp_eval (u, space, geometry, vtk_pts);
% [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
% figure ('Units', 'pixels', 'Position', [150 200 1000 350])
% subplot (1,2,1)
% h1 = pcolor (X, Y, real(eu));
% colorbar
% h1.EdgeColor = 'none';
% h1.FaceColor = 'interp';
% title ('Numerical solution: \Re(u_h)'), axis tight
% xlabel('x')
% ylabel('Time')
% subplot (1,2,2)
% h2 = pcolor (X, Y, real(problem_data.uex (X,Y)));
% colorbar
% h2.EdgeColor = 'none';
% h2.FaceColor = 'interp';
% title ('Exact solution: \Re(u)'), axis tight
% xlabel('x')
% ylabel('Time')
% 
%% 5) DISPLAY ERRORS of the computed solution in the L2 and H1 norm
% compute the error for the real part:
Uex = @(x, y, t) real(problem_data.uex(x, y, T*t));
GradUex = @(x, y, t) real(problem_data.graduex(x, y, T*t));
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, real(u), Uex, GradUex)


%% 6) Save solution
n = method_data.nsub(1);
d = method_data.degree(1);
filename = ['test_schrodinger_degree_' num2str(d) '_subs_' num2str(n) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);

