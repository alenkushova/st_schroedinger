%% WITH L^2 NORM
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 3
    err_l2 = [];  NN = [];
for i = 1 : 7
    N = 2^(i+2);
    % load the workspace
    s = ['test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2_new];

end
    deg = ['p = ' num2str(grad) ', uniform ref.'];
    loglog(NN,err_l2,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
%    loglog(NN,10^3*NN.^(grad+1),'-.k','Linewidth',1.5,'DisplayName','h^{4}')
    clear d s
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(deg(1)+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('||u-u_h||_{L^2}')

badmeshes=[5 17 42 70 110];

for grad = 3
    err_l2 = [];  HH = [];
for i = 1 : 5
    N = badmeshes(i);
    % load the workspace
    s = ['test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    HH = [HH 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2_new];
end
    deg = ['p = ' num2str(grad) ', bad meshes'];
    loglog(HH,err_l2,'-*','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    loglog(NN,10^3*NN.^(4),'-.k','Linewidth',1.5,'DisplayName','h^{2}')
    hold off
    clear d s
end

%% WITH GRAPH NORM
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 3
    err_G = [];  NN = [];
for i = 1 : 7
    N = 2^(i+2);
    % load the workspace
    s = ['Test_Graph_convergence/test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_G  = [err_G error_Graph];
end
    deg = ['p = ' num2str(grad) ', uniform ref.'];
    loglog(NN,err_G,'-o','Linewidth',1.5,'DisplayName',deg)
%    grid on, hold on
%    loglog(NN,10^5*NN.^(2),'-.k','Linewidth',1.5,'DisplayName','h^{2}')
    clear d s
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(deg+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('||u-u_h||_{Graph}')

badmeshes=[5 17 42 70 110];

for grad = 3
    err_G = [];  HH = [];
for i = 1 : 5
    N = badmeshes(i);
    % load the workspace
    s = ['Test_Graph_conv_badmeshes/test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    HH = [HH 1/N];
%         % save errors: this saves the error i already computed
    err_G  = [err_G error_Graph];
end
    deg = ['p = ' num2str(grad) ', bad meshes'];
    loglog(HH,err_G,'-*','Linewidth',1.5,'DisplayName',deg)
%    grid on, hold on
    loglog(NN,10^5*NN.^(2),'-.k','Linewidth',1.5,'DisplayName','h^{2}')
    hold off
    clear d s
end
% legend('Location','southeast')
% loglog(NN,(6*1./sqrt(NN)).^(deg+1) ,'-s','Linewidth',1.5)
% title('Error convergence for solution with Fast Diagonalization')
% xlabel('h')
% ylabel('||u-u_h||_{Graph}')
