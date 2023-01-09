%% Test with degree 3 splines and G - l2 norm
p = 3; eigval = zeros(1,5); infsup = zeros(1,5); h = zeros(1,5);
geometries = {}; meshes = {}; spaces = {}; eigvect = {}; lamMww = {};
v = {}; Av = {};
% 'filename' of the file where to save the results
filename = ['degree_' num2str(p) '_infsup_Gsn_l2_results.mat'];
% load(filename)
for i = 1:128
  nel = i+4; j = i;
%  [mu, ei, geo, msh, space, w, tmp_lamMww, tmpv, tmpAv]  = infsuptest (p, nel, 'l2', 'l2');
  [mu, ei, geo, msh, space, w]  = infsuptest (p, nel, 'G', 'l2');
  h (j) = 1/nel;   infsup (j) = mu;   eigval (j) = ei;
  eigvect    = cat(2, eigvect, {w});
%  lamMww     = cat(2, lamMww, {tmp_lamMww});
%  v          = cat(2, v, {tmpv});
%  Av         = cat(2, Av, {tmpAv});
  geometries = cat(2, geometries, {geo});
  meshes     = cat(2, meshes, {msh});
  spaces     = cat(2, spaces, {space});
% Saving at any loop iteration:
  save( filename , 'h', 'infsup', 'eigval', 'eigvect', 'geometries', 'meshes', 'spaces')
  fprintf(['Iteration ' num2str(i) '. Number of elements ' num2str(nel) '. Results saved in file: ' filename '\n\n'])
end
