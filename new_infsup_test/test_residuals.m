%% Test with degree 3 splines
p = 3; 
h = zeros(1,5);
geometries = {};
meshes     = {};
spaces     = {};
residuals  = {};

for i = 1:100
  nel = i+4;
  j = i;
  [geo, msh, sp, res] = residual_rates(p, nel);
%  [mu, ei, geo, msh, space, w, tmp_lamMww, tmpv, tmpAv]  = infsuptest (p, nel);
%  [mu, ei, geo, msh, space, w, tmpv]  = infsuptest (p, nel);
  geometries = cat(2, geometries, {geo});
  meshes     = cat(2, meshes, {msh});
  spaces     = cat(2, spaces, {sp});
  residuals  = cat(2, residuals, {res});
  h (j) = 1/nel;
end

%% plot residuals:
loglog(h,[residuals{:}],'LineWidth',1.5)
grid on
xlabel('h')
legend('|| (Id-PW) Sv ||_{L^2} / ||w||_{L^2} ','location','northeast')
title('residuals sequence with respect to mesh size')

%% save
filename = ['residuals_degree_' num2str(p) '.mat'];
save( filename , 'h', 'geometries', 'meshes', 'spaces', 'residuals')
fprintf (['Results saved in file: ' filename '\n\n'])
