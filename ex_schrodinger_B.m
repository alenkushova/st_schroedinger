%% EX_SCHRODINGER_B: solve the Schroedinger problem in one dimension. (S-T)
clear 
close all
clc

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data  
% Physical domain, defined as NURBS map given in a text file
problem_data.geo_name = nrbline ([0 0], [1 0]);

% Type of boundary conditions for each side of the domain

problem_data.nmnn_sides   = []; % Neumann 
problem_data.drchlt_sides = [1 2];  % Dirichlet
problem_data.init_sides   = [1];  % initial value sides (for time direction)
problem_data.prdc_sides   = []; % Periodic


% Parameters:
omg = 0.2;
a0  = (2/(omg^2))^(1/4);

% Exact solution:
problem_data.uex     = @(x, t) (a0*exp(-1i*(x.^2+t.^2)/(omg^2)));
problem_data.graduex = @(x, t) (cat (1, ...
                reshape (-2i*x/(omg^2).*problem_data.uex(x,t) , [1, size(x)]), ...
                reshape (-2i*t/(omg^2).*problem_data.uex(x,t) , [1, size(x)])));

% Source term
problem_data.f = @(x, t) ...
            (2i*omg^2 + 4*x.^2 + 2*omg^2*t)/omg^4.*problem_data.uex(x,t);
problem_data.h = @(x, t, ind) problem_data.uex(x,t);

%in this moment this values are fixed! Need to define better solver
problem_data.dimension  = 1;      % Space Dimension
problem_data.T          = 1;      % Final time

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = 3;        % Degree of the splines
method_data.regularity = method_data.degree-1; % Regularity of the splines
method_data.nsub       = 32;        % Number of subdivisions
method_data.nquad      = method_data.degree+1; % Points for the Gaussian quadrature rule
method_data.solver     = 'FD';      % Solver: Fast Diag 'FD' or Matlab Backslash 'M'

%% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_schrodinger_st (problem_data, method_data);

%% 4) POSTPROCESSING
vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)};
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
subplot (1,2,1)
h1 = pcolor (X, Y, real(eu));
colorbar
h1.EdgeColor = 'none';
h1.FaceColor = 'interp';
title ('Numerical solution: \Re(u_h)'), axis tight
xlabel('x')
ylabel('Time')
subplot (1,2,2)
h2 = pcolor (X, Y, real(problem_data.uex (X,Y)));
colorbar
h2.EdgeColor = 'none';
h2.FaceColor = 'interp';
title ('Exact solution: \Re(u)'), axis tight
xlabel('x')
ylabel('Time')

%% 5) DISPLAY ERRORS of the computed solution in the L2 and H1 norm
% compute the error for the real part:
Uex = @(x,t) real(problem_data.uex(x,t));
GradUex = @(x, t) real(problem_data.graduex(x,t));
[error_h1, error_l2] = ...
           sp_h1_error (space, msh, real(u), Uex, GradUex)


%% 6) Save solution
n = method_data.nsub;
d = method_data.degree;
filename = ['test_schrodinger_degree_' num2str(d) '_subs_' num2str(n) '.mat'];
save(filename)
fprintf ('The result is saved in the file: %s \n \n', filename);

