%% Analyze the errors into a plot!
close all
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for d = 1 : 7
    err_l2 = [];  NN = [];
for i = 1 : 5
    N = 2^(i+2);
    % load the workspace
    s = ['Fast Diagonalization with Schur\test_schrodinger_degree_' num2str(d) '_subs_' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2];
end
    deg = ['p = ' num2str(d)];
    loglog(NN,err_l2,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(d+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Fast Diagonalization')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|')
%% Second figure
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
for d = 1 : 7
    err_l2 = [];  NN = [];
for i = 1 : 5
    N = 2^(i+2);
    % load the workspace
    s = ['Matlab backslash\test_schrodinger_degree_' num2str(d) '_subs_' num2str(N) '.mat'];    
    load(s);
    NN = [NN 1/N];
%         % save errors: this saves the error i already computed
    err_l2  = [err_l2 error_l2];
end
    deg = ['p = ' num2str(d)];
    loglog(NN,err_l2,'-o','Linewidth',1.5,'DisplayName',deg)
    grid on, hold on
end
legend('Location','southeast')
%loglog(NN,(6*1./sqrt(NN)).^(d+1) ,'-s','Linewidth',1.5)
title('Error convergence for solution with Matlab backslash')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|')