function errl2 = schroedinger_l2_error_projection(spacev, msh, u, spacew, Psu)
  errl2 = 0;
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    spv_col  = sp_evaluate_col (spacev, msh_col, 'value', true,...
        'gradient', true, 'hessian', true);
    spw_col = sp_evaluate_col (spacew, msh_col, 'value', true, 'gradient', false);
    errl2 = errl2 + (schroedinger_l2_error (spv_col, msh_col, u, spw_col, Psu)).^2;
  end
  
  errl2 = sqrt (errl2);

end

function [err, err_elem] = schroedinger_l2_error (spv, msh, u, spw, Psu)

  valgradu = sp_eval_msh (u, spv, msh, 'gradient');
  valdtu = reshape(valgradu(end,:,:), spv.ncomp, msh.nqn, msh.nel);
    
  valhessu = sp_eval_msh (u, spv, msh, 'hessian');
  vallapu = zeros(spv.ncomp, msh.nqn, msh.nel);
  for ii = 1: msh.rdim-1
    vallapu = vallapu + reshape(valhessu(ii,ii,:,:),spv.ncomp, msh.nqn, msh.nel);
  end
    
  valsu  = 1i*valdtu - vallapu;
   
  valPsu = sp_eval_msh (Psu, spw, msh);
  valPsu = reshape(valPsu, spv.ncomp, msh.nqn, msh.nel);
       
  w = msh.quad_weights .* msh.jacdet;
    
  errl2_elem = sum(reshape(sum((abs(valsu - valPsu)).^2, 1), [msh.nqn, msh.nel]) .* w);
  
  err  = sqrt (sum (errl2_elem));
  err_elem  = sqrt (errl2_elem);

end