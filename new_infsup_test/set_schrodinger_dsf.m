% this function sets the discrete tp B-spline spaces for the schroedinger
% problem, that is S^{p+1,p+2} space-time domain of maximum smoothness both for
% trial and S^{p,p} for test functions.
%
% Call: [arguments] = set_schrodinger_dwf (trial_p, test_p, nel)
%
% Input:
%        nel = number of elements per univariate direction
%    trial_p = trial polinomial degree
%     test_p = test polinomial degree
%
% Output: 'output_args' structure with the following fields:
%    geometry = geometry of the domain
%         msh = mesh associated to B-splines space-time 
%        mshx = mesh associated to B-splines in x direction(s)
%        msht = mesh associated to B-splines in time direction 
%       space = B-spline space 
%     x_space =   "        "   in x direction(s)
%     t_space =   "        "   in time direction
%

function output_args = set_schrodinger_dsf (nel, trial_p, test_p)
T = 1; 
% Geometries.
output_args.xt_geo_name = 'geo_square.txt';
output_args.x_geo_name  = nrbline ([0 0], [1 0]);
output_args.t_geo_name  = nrbline ([0 0], [T 0]);

% Choice of discretization parameters. 
output_args.trial_degree     = [trial_p+2 trial_p+1];                                % Degree of the splines (last is time dir)
output_args.trial_regularity = output_args.trial_degree-1;                             % Regularity of the splines ( " )
output_args.test_degree      = [test_p+2 test_p+1];                                % Degree of the splines (last is time dir)
output_args.test_regularity  = [test_p-1 test_p-1];                             % Regularity of the splines ( " )
output_args.nsub  = [nel T*nel];                                      % Number of subdivisions ( " )
output_args.nquad = max(output_args.trial_degree, output_args.test_degree) +1;                             % Points for the Gaussian quadrature rule ( " )
end