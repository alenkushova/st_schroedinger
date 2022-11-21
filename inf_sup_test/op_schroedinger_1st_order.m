% This function assambles the matrix associated to the schroedinger 
% operator exploiting the tensor-product structure of B-spline spaces.
% 
% Call:               
%                A = op_schroedinger_1st_order (mshx, msht, test_spx, test_spt, options)
% [rows,cols,vals] = op_schroedinger_1st_order (mshx, msht, test_spx, test_spt, options)
% 
% Input:       
%              mshx = mesh for x direction(s)
%              msht = mesh for time direction
%         trial_spx = trial space in x-direction(s)
%         trial_spt = tiral space in time direction
%           options = you can add a different space for test functions 
%                
%     I.E. test_spx = spazio delle funzioni test in x diverso dalle trial
%          test_spx = spazio delle funzioni test in tempo diverso dalle trial
%
% Output:    A = matrix of dicrete Schroedinger operator associated to the mesh.
% 
function varargout = op_schroedinger_1st_order (mshx, msht, trial_spx, trial_spt,...
                                           test_spx, test_spt)
      
  if (nargin < 4)
    error('op_schroedinger_1st_order: not enought imput arguments');
  elseif (nargin == 4)
    test_spx = trial_spx;
    test_spt = trial_spt;
  elseif (nargin == 5 ) 
    test_spt = trial_spt;
  end
      
  Wt = op_gradu_v_tp (trial_spt, test_spt, msht);
  Mt = op_u_v_tp (trial_spt, test_spt, msht);
  Ms = op_u_v_tp (trial_spx, test_spx, mshx);
  Ks = op_gradu_gradv_tp (trial_spx, test_spx, mshx);
  A  = 1i*kron(Wt',Ms) + kron(Mt,Ks);
  
  if (nargout == 1 || nargout == 0)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, values] = find(A);
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_schroedinger_1st_order: wrong number of output arguments')
  end

end


