%% Analyze the errors into a plot! L2-convergence for fixed M = 1250
% Space-time refinement with k = sqrt(h) and h = space meshwidth, k = time
% meshwidth.

names = {'Eulero esplicito','Crank Nicolson','Eulero implicito'}; 
rates = {@(h) 2.2*h.^(0.25),@(h) 1.7*h.^(0.33),@(h) 2.2*h.^(0.25)};
rate_names = {'h^{0.25}','h^{0.5}','h^{0.25}'};
for name = 2
    theta = (name-1)/2;
    conv_rate = rates{name};
    errore_l2 = [];  NN = [];
for i = 1 : 5
    subdivisions = 2^(i+1);
    % load the workspace
    s = [names{name} '\test_schrodinger_degree_3_subs_'...
           num2str(subdivisions) '_theta_' num2str(theta) '.mat'];    
    load(s);
    NN = [NN 1/subdivisions];
%         % save errors: this saves the error i already computed
    errore_l2  = [errore_l2 errore_l2_L2];

end
    figure ('Units', 'pixels', 'Position', [150 200 1000 350])
    loglog(NN,errore_l2,'-o','Linewidth',1.5,'DisplayName',names{name})
    grid on, hold on
    legend('Location','southeast')
    loglog(NN,conv_rate(NN) ,'-sk','Linewidth',1.5,'DisplayName',rate_names{name} )
    title('Error convergence for solution with Theta-method')
    xlabel('h')
    ylabel('|\Re(u)-\Re(u_h)|_{L^2}')
    clear s
end


%% Analyze the errors into a plot! L2-convergence for fixed M = 1250
% Time refinement with h = space meshwidth, h = 1/128, k = time meshwidth.

names = {'Crank Nicolson','Eulero implicito'}; 
rates = {@(h) 2.7*h.^(0.25),@(h) 3.2*h.^(0.25)};
rate_names = {'h^{0.25}','h^{0.25}'};
for name = 1 : 2
    theta = name/2;
    conv_rate = rates{name};
    errore_l2 = [];  NN = [];
for i = 1 : 5
    subdivisions = 2^(i+1);
    % load the workspace
    s = [names{name} '\test_schrodinger_degree_3_subs_128_and_'...
           num2str(2*subdivisions) '_theta_' num2str(theta) '.mat'];    
    load(s);
    NN = [NN 1/subdivisions];
%         % save errors: this saves the error i already computed
    errore_l2  = [errore_l2 errore_l2_L2];

end
    figure ('Units', 'pixels', 'Position', [150 200 1000 350])
    loglog(NN,errore_l2,'-o','Linewidth',1.5,'DisplayName',names{name})
    grid on, hold on
    legend('Location','southeast')
    loglog(NN,conv_rate(NN) ,'-sk','Linewidth',1.5,'DisplayName',rate_names{name} )
    title('Error convergence for solution with Theta-method')
    xlabel('h')
    ylabel('|\Re(u)-\Re(u_h)|_{L^2}')
    clear s
end
