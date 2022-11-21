%% Test with degree 3 splines
p = 3; 
eigval = zeros(1,5);
eigvect = [];
infsup = zeros(1,5);
h = zeros(1,5);
geometries = {};
meshes = {};
spaces = {};

for i = 1:2
  nel = 10*i;
  j = i;
  [mu, ei, geo, msh, space, w]  = infsuptest (p, nel);
  infsup (j) = mu;
  eigval (j) = ei;
  eigvect = [eigvect w]; 
  geometries = cat(2, geometries, {geo});
  meshes     = cat(2, meshes, {msh});
  spaces     = cat(2, spaces, {space});
  h (j) = 1/nel;
end
%% plot infsup:
loglog(h,infsup,'LineWidth',1.5)
grid on
hold on
loglog(h,h.^(1/3),':','LineWidth',1.5)
loglog(h,h.^(2/3),'-.','LineWidth',1.5)
xlabel('h')
legend('\mu_{min}^{1/2}','h^{1/3}','h^{2/3}','location','northwest')
title('Inf-sup test')
hold off

%% plot eigenvectors:

%% save
filename = ['new_degree_' num2str(p) '_infsup_graphnorm_results.mat'];
save( filename , 'h', 'infsup', 'eigenv')
fprintf (['Results saved in file: ' filename '\n\n'])