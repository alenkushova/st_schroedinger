% post-processing of the inf-sup conditions.
clear 
close all
clc

load('database/Nuovo_degree_3_infsup_l2_l2_results.mat')
loglog(h,infsup,'LineWidth',1.5)
grid on
hold on
xlabel('h')
title('Inf-sup test')

load('database/Nuovo_degree_3_infsup_l2_G_results.mat')
loglog(h,infsup,'LineWidth',1.5)

load('database/Nuovo_degree_3_infsup_G_l2_results.mat')
loglog(h,infsup,'LineWidth',1.5)

loglog(h,0.6*h.^(1/2),'--k','LineWidth',1.5)
legend('L2-L2','Graph-L2','L2-Graph','h^{1/2}','location','southeast')

hold off

