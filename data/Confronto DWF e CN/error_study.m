%% Analyze the errors into a plot! L2-convergence for fixed M = 1250
% Space-time refinement with k = sqrt(h) and h = space meshwidth, k = time
% meshwidth.

theta = 1/2; folder = 'spazio_fine';                                       % Crank Nicolson + spazio fine
% theta = 1/2; folder = 'tempo_fine';                                        % Crank Nicolson + tempo fine
names = {'TM','DWF','L2'};
TM_u_errore_l2 = []; ST_u_errore_l2 = []; L2_u_errore_l2 = [];  NN = [];
for i = 4 : 8
    subs = 2^(i);
    % load the workspace
    s = [ folder '\test_schrodinger_degree_3_subs_'...
           num2str(subs) '_and_' num2str(2*subs) '_theta_' num2str(theta) '.mat'];    
    load(s);
    NN = [NN 1/subs];
%         % save errors: this saves the error i already computed
    TM_u_errore_l2  = [TM_u_errore_l2 errore_l2_L2];
    ST_u_errore_l2  = [ST_u_errore_l2 u_error_l2];
    L2_u_errore_l2  = [L2_u_errore_l2 Pu_error_l2];
    
end

figure ('Units', 'pixels', 'Position', [150 200 1000 350])
loglog(NN,TM_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{1})
grid on, hold on
loglog(NN,ST_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{2})
loglog(NN,L2_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{3})
legend('Location','southeast')
loglog(NN,2.2*NN.^(0.25),'-sk','Linewidth',1.5,'DisplayName','h^{0.25}' )
loglog(NN,2.6*NN.^(0.2),'-.sk','Linewidth',1.5,'DisplayName','h^{0.2}' )
title('Crank Nicolson, DWF and L2-projection with h = dt')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|_{L^2}')
axis tight


%% draw slopes on boxes
TM_slopes = log(TM_u_errore_l2(1:end-1)./TM_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
ST_slopes = log(ST_u_errore_l2(1:end-1)./ST_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
L2_slopes = log(L2_u_errore_l2(1:end-1)./L2_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
x_points  = (NN(1:end-1)+NN(2:end))/2;
TM_points = (TM_u_errore_l2(1:end-1)+TM_u_errore_l2(2:end))/2;
ST_points = (ST_u_errore_l2(1:end-1)+ST_u_errore_l2(2:end))/2;
L2_points = (L2_u_errore_l2(1:end-1)+L2_u_errore_l2(2:end))/2;
for i = 4
xpnt= x_points(i); ypnt = TM_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .01/(2^(i-1)) .13];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.8)
text(xpnt,ypnt,num2str(TM_slopes(i)));
end

for i = 4
xpnt= x_points(i); ypnt = ST_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .01/(2^(i-1)) .08];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.9/1.2^(i-1))
text(xpnt,ypnt,num2str(ST_slopes(i)));
end

for i = 4
xpnt= x_points(i); ypnt = L2_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .01/(2^(i-1)) .09/1.2^(i-1)];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.8)
text(xpnt,ypnt,num2str(L2_slopes(i)));
end

%% Analyze the errors into a plot! L2-convergence for fixed M = 1250
% Space-time refinement with k = sqrt(h) and h = space meshwidth, k = time
% meshwidth.

theta = 1/2; folder = 'tempo_fine';                                         % Eulero implicito + spazio fine
% theta = 1; folder = 'tempo_fine';                                         % Eulero implicito + tempo fine
names = {'TM','DWF','L2'};
TM_u_errore_l2 = []; ST_u_errore_l2 = []; L2_u_errore_l2 = [];  NN = [];
for i = 3 : 7
    subs = 2^(i);
    % load the workspace
    s = [ folder '\test_schrodinger_degree_3_subs_'...
           num2str(subs) '_and_1024_theta_' num2str(theta) '.mat'];    
    load(s);
    NN = [NN 1/subs];
%         % save errors: this saves the error i already computed
    TM_u_errore_l2  = [TM_u_errore_l2 errore_l2_L2];
    ST_u_errore_l2  = [ST_u_errore_l2 u_error_l2];
    L2_u_errore_l2  = [L2_u_errore_l2 Pu_error_l2];
    
end
%%
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
loglog(NN,TM_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{1})
grid on, hold on
loglog(NN,ST_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{2})
loglog(NN,L2_u_errore_l2,'-o','Linewidth',1.5,'DisplayName',names{3})
legend('Location','southeast')
loglog(NN,15.e-6*NN.^(0.5),'-sk','Linewidth',1.5,'DisplayName','h^{0.5}' )
title('Crank Nicolson, DWF and L2-projection with dt << h^2')
xlabel('h')
ylabel('|\Re(u)-\Re(u_h)|_{L^2}')

%% draw slopes on boxes
TM_slopes = log(TM_u_errore_l2(1:end-1)./TM_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
ST_slopes = log(ST_u_errore_l2(1:end-1)./ST_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
L2_slopes = log(L2_u_errore_l2(1:end-1)./L2_u_errore_l2(2:end))./log(NN(1:end-1)./NN(2:end));
x_points  = NN(2:end);
TM_points = TM_u_errore_l2(2:end);
ST_points = ST_u_errore_l2(2:end);
L2_points = L2_u_errore_l2(2:end);

formatSpec = '%.2f';

for i = 4
xpnt= x_points(i); ypnt = TM_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .001 1.3e-7];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.8)
text(xpnt,ypnt,num2str(TM_slopes(i), formatSpec));
end

for i = 4
xpnt= x_points(i); ypnt = ST_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .001 1.3e-7];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.9/1.2^(i-1))
text(xpnt,ypnt,num2str(ST_slopes(i), formatSpec));
end

for i = 4
xpnt= x_points(i); ypnt = L2_points(i);
pos = [(xpnt/1.015) (ypnt/1.045) .001 1.3e-7];
rectangle('Position',pos,'FaceColor','#FFFF00','LineWidth',0.8)
text(xpnt,ypnt,num2str(L2_slopes(i), formatSpec));
end
