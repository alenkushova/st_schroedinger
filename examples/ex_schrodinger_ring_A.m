%% EX_SCHRODINGER: solve the Schroedinger problem in one dimension. (S-T)
clear 
close all
clc
% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
T = 1; 
problem_data.T = T ;      % Final time. The geometry supports just T=1 for the moment. 

% Physical domain, defined as NURBS map given in a text file
problem_data.xt_geo_name = 'geo_ring_time.txt'; % goemetry file of the space-time domain! 
problem_data.x_geo_name = 'geo_ring.txt';
% problem_data.xt_geo_name = 'geo_cube.txt'; % goemetry file of the space-time domain! 
% problem_data.x_geo_name = 'geo_square.txt';
problem_data.t_geo_name = nrbline ([0 0], [T 0]);

% Type of boundary conditions for each side of the domain
problem_data.nmnn_sides     = []; % Neumann 
problem_data.drchlt_sides   = [1 2 3 4 5];  % Dirichlet
problem_data.x_drchlt_sides = [1 2 3 4];  % Dirichlet
problem_data.prdc_sides     = []; % Periodic

% Parameters:
Mm  = 1.5;
T_0 = 1.5;
beta= 2.5;
R = @(t) (T_0^2 - 1i*beta*t);
DR= -1i*beta;

% Exact solution:
problem_data.uex     = @(x, y, t) (1./R(t).*exp(-(x.^2 + y.^2)./R(t)));
problem_data.graduex = @(x, y, t) (cat (1, ...
                reshape (-(2*x./R(t)).*problem_data.uex(x, y, t) , [1, size(x)]), ...
                reshape (-(2*y./R(t)).*problem_data.uex(x, y, t) , [1, size(y)]), ...
                reshape (((x.^2+y.^2)./R(t) - 1).*(DR./R(t)).*problem_data.uex(x, y, t) , [1, size(t)])));
% Source term
problem_data.f = @(x, y, t) ...
            ((x.^2+y.^2)./R(t) - 1).*((-4+1i*DR)./R(t)).*problem_data.uex(x, y, t);
problem_data.h = @(x, y, t, ind) problem_data.uex(x, y, t);
problem_data.x_h = @(x, y, ind) zeros (size (x)); % auxiliary function. 
problem_data.gmm = 1i;
problem_data.eta = 1;


% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3 3]; % Degree of the splines (last is time dir)
method_data.regularity = method_data.degree-1; % Regularity of the splines
method_data.nsub       = [64 64 64]; % Number of subdivisions
method_data.nquad      = method_data.degree+1; % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';     % Fast Diag 'FD' or Matlab Backslash 'M'

% 3) CALL TO THE SOLVER
[geometry, msh, space, u, gmres_output] = solve_schrodinger_gmres (problem_data, method_data);

% 4) POST-PROCESSING
% 4.1) EXPORT TO PARAVIEW
% output_file = 'Schroedinger_st_2D';
% vtk_pts = {linspace(0, 1, 200), linspace(0, 1, 200), linspace(0, 1, 20)};
% fprintf ('The result is saved in the file %s \n \n', output_file);
% sp_to_vtk (u, space, geometry, vtk_pts, output_file, 'u')

% 4.2) FIGURE IN MATLAB
figure
sp_plot_solution (real(u), space, geometry, [40 40 40], [20 20 20]);
shading interp 

gmres_output

% 5) DISPLAY ERRORS of the computed solution in the L2 and GRAPH norm
Uex = @(x, y, t) (problem_data.uex(x, y, t));
rhs = @(x, y, t) (problem_data.f(x, y, t));
[error_Graph, error_l2_new] = schroedinger_graph_error(space, msh, u, Uex, rhs)

% 6) Save solution
n = method_data.nsub(1);
d = method_data.degree(1);
filename = ['test_schrodinger_gmres_degree_' num2str(d) '_subs_' num2str(n) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);

