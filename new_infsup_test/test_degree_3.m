%% Test with degree 3 splines
p = 3; 
eigval = zeros(1,5);
infsup = zeros(1,5);
h = zeros(1,5);
geometries = {};
meshes = {};
spaces = {};
eigvect = {};
lamMww = {};
v = {};
Av = {};

for i = 1:100
  nel = i+4;
  j = i;
  [mu, ei, geo, msh, space, w, tmp_lamMww, tmpv, tmpAv]  = infsuptest (p, nel);
%  [mu, ei, geo, msh, space, w, tmpv]  = infsuptest (p, nel);
  infsup (j) = mu;
  eigval (j) = ei;
  eigvect    = cat(2, eigvect, {w});
  lamMww     = cat(2, lamMww, {tmp_lamMww});
  v          = cat(2, v, {tmpv});
  Av         = cat(2, Av, {tmpAv});
  geometries = cat(2, geometries, {geo});
  meshes     = cat(2, meshes, {msh});
  spaces     = cat(2, spaces, {space});
  h (j) = 1/nel;
end

%% plot infsup:
loglog(h,infsup,'LineWidth',1.5)
grid on
hold on
loglog(h,0.6*h.^(1/3),':','LineWidth',1.5)
loglog(h,0.6*h.^(2/3),'-.','LineWidth',1.5)
xlabel('h')
legend('\mu_{min}^{1/2}','h^{1/3}','h^{2/3}','location','northwest')
title('Inf-sup test')
hold off

%% plot eigenvectors:
vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)};
for i = 1 
    [ew, ~] = sp_eval (eigvect{i}, spaces{i}, geometries{i}, vtk_pts);
    [elamMww, ~] = sp_eval (lamMww{i}, spaces{i}, geometries{i}, vtk_pts);
    [ev, ~] = sp_eval (v{i}, spaces{i}, geometries{i}, vtk_pts);
    [eAv, F] = sp_eval_schroedinger (v{i}, spaces{i}, geometries{i}, vtk_pts, 'schroedinger');
%    [eAv, F]= sp_eval (Av{i}, spaces{i}, geometries{i}, vtk_pts);
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    figure ('Units', 'pixels', 'Position', [50 200 1200 350])
    subplot (1,3,1)
    h1 = pcolor (X, Y, real(ew));
    colorbar
    colormap jet
    h1.EdgeColor = 'none';
    h1.FaceColor = 'interp';
    title ('Eigenvector w_{min}'), axis tight
    xlabel('x')
    ylabel('Time')
    subplot (1,3,2)
    h2 = pcolor (X, Y, real(ev));
    colorbar
    colormap jet
    h2.EdgeColor = 'none';
    h2.FaceColor = 'interp';
    title ('v_{min} = (Mv^{-1})A^{*}w_{min}'), axis tight
    xlabel('x')
    ylabel('Time')
    subplot (1,3,3)
    h3 = pcolor (X, Y, real(eAv-elamMww));
    colorbar
    colormap jet
    h3.EdgeColor = 'none';
    h3.FaceColor = 'interp';
    title (' Av_{min} - \lambda_{min}^2 Mw w_{min}'), axis tight
    xlabel('x')
    ylabel('Time')
end


%% save
filename = ['degree_' num2str(p) '_infsup_l2_l2_results.mat'];
save( filename , 'h', 'infsup', 'eigval', 'eigvect', 'v', 'Av', 'geometries', 'meshes', 'spaces')
fprintf (['Results saved in file: ' filename '\n\n'])