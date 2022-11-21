% post-processing of the inf-sup conditions.
clear 
close all
clc

load('degree_3_infsup_l2_l2_results.mat')
loglog(h,infsup,'LineWidth',1.5)
grid on
hold on
xlabel('h')
title('Inf-sup test')

load('degree_3_infsup_l2_graphnorm_results.mat')
loglog(h,infsup,'LineWidth',1.5)

load('degree_3_infsup_graphnorm_l2_results.mat')
loglog(h,infsup,'LineWidth',1.5)
legend('L2-L2','Graph-L2','L2-Graph','location','southeast')
hold off