%% Analyze the errors into a plot! L2-convergence for fixed M = 625
close all
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 1 : 5
    err_l2 = [];  NN = [];
for i = 1 : 7
    N = 2^(i+2);
    % load the workspace
    s = ['Truncated expansion for M625 lq\test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(2*N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2];

end
    deg = ['p = ' num2str(grad)];
    loglog(NN,err_l2,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    clear grad s
end
legend('Location','southeast')
loglog(NN,1.7*NN.^(0.25) ,'-sk','Linewidth',1.5,'DisplayName','h^{0.2}')
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|_{L^2}')


%% H1-convergence for fixed M = 625
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for grad = 1 : 5
    err_h1 = [];  NN = [];
for i = 1 : 7
    N = 2^(i+2);
    % load the workspace
    s = ['Truncated expansion for M625 lq\test_schrodinger_degree_' num2str(grad)...
           '  ' num2str(grad) '_subs_' num2str(N) '  ' num2str(2*N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_h1  = [err_h1 error_h1];
end
    deg = ['p = ' num2str(grad)];
    loglog(NN,err_h1,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
    clear d s
end
legend('Location','southwest')
%loglog(NN,(6*1./sqrt(NN)).^(d+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|_{H^1}')
