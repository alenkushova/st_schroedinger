%% Analyze the errors into a plot! L2-convergence 
badmeshes=[5 17 42 70 110];
close all
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 3
    err_l2 = [];  NN = [];
for i = 1 : 5
    N = badmeshes(i);
    % load the workspace
    s = ['test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2_new];

end
    deg = ['p = ' num2str(grad)];
    loglog(NN,err_l2,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    loglog(NN,10^3*NN.^(grad+1),'-.k','Linewidth',1.5,'DisplayName','h^{4}')
    hold off
    clear d s
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(deg(1)+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('||u-u_h||_{L^2}')


%% Graph norm-convergence 
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 3
    err_G = [];  NN = [];
for i = 1 : 5
    N = badmeshes(i);
    % load the workspace
    s = ['test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_G  = [err_G error_Graph];
end
    deg = ['p = ' num2str(grad)];
    loglog(NN,err_G,'-*','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    loglog(NN,10^5*NN.^(2),'-.k','Linewidth',1.5,'DisplayName','h^{2}')
    hold off
    clear d s
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(deg+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('||u-u_h||_{Graph}')