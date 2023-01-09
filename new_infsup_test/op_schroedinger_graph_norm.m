% This function assambles the matrix associated to the schroedinger 
% operator graph norm.
% 
% Call:  Mv = op_schroedinger_graph_norm (msh, mshx, msht, space, spx, spt)
% 
% Input:       
%         msh   = mesh for space-time grid
%         mshx  = mesh for x direction(s)
%         msht  = mesh for time direction
%         space = space-time B-splines
%         spx   = space in x direction(s)
%         spt   = space in time direction
%                
%     I.E. test_spx = spazio delle funzioni test in x diverso dalle trial
%          test_spx = spazio delle funzioni test in tempo diverso dalle trial
%
% Output: Mv = matrix representing the discrete graph norm associated to
%              Schroedinger operator in space-time discretizations.
% 
function Mv = op_schroedinger_graph_norm (msh, mshx, msht, space, spx, spt)
        
  Mt     = op_u_v_tp (spt, spt, msht);
  Ms     = op_u_v_tp (spx, spx, mshx);
  dtdt   = op_gradu_gradv_tp (spt, spt, msht);
  laplap = op_laplaceu_laplacev_tp (spx, spx, mshx);
  dtlap  = op_dtu_laplacev_tp (space, space, msh);
  idtlap = -1i*dtlap;
  mat = kron(dtdt, Ms) + idtlap + idtlap' + kron(Mt,laplap);
%  Mv  = op_u_v_tp (space, space, msh) + mat; % Graph norm
  Mv  = mat;  % Graph seminorm
 
end


