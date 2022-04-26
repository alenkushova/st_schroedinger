clear
close all
Mod = 5;
p1= 3;
p2= 4;
h1= 128;
h2= 256;
err_l2_deg_3_128 = zeros(1,4); err_h1_deg_3_128 = zeros(1,4);
err_l2_deg_3_256 = zeros(1,4); err_h1_deg_3_256 = zeros(1,4);
err_l2_deg_4_128 = zeros(1,4); err_h1_deg_4_128 = zeros(1,4);
err_l2_deg_4_256 = zeros(1,4); err_h1_deg_4_256 = zeros(1,4);
modes = [5 25 125 625];
for i = 1 : 4
%__________________________________________________________________________    
    str1 = ['Truncated expansion for M' num2str(Mod^i)...
        '\test_schrodinger_degree_3  3_subs_128  256.mat'];
    load(str1);
    err_l2_deg_3_128(i) = error_l2; err_h1_deg_3_128(i) = error_h1;
%__________________________________________________________________________    
    str2 = ['Truncated expansion for M' num2str(Mod^i)...
        '\test_schrodinger_degree_3  3_subs_256  512.mat'];
    load(str2);
    err_l2_deg_3_256(i) = error_l2; err_h1_deg_3_256(i) = error_h1;
%__________________________________________________________________________    
    str3 = ['Truncated expansion for M' num2str(Mod^i)...
        '\test_schrodinger_degree_4  4_subs_128  256.mat'];
    load(str3);
    err_l2_deg_4_128(i) = error_l2; err_h1_deg_4_128(i) = error_h1;
%__________________________________________________________________________    
    str4 = ['Truncated expansion for M' num2str(Mod^i)...
        '\test_schrodinger_degree_4  4_subs_256  512.mat'];
    load(str4);
    err_l2_deg_4_256(i) = error_l2; err_h1_deg_4_256(i) = error_h1;
end
figure ('Units', 'pixels', 'Position', [150 200 1000 350])
hold on; grid on;
loglog(modes,err_l2_deg_3_128,'-o','Linewidth',1.5,'DisplayName','p = 3, subs = 128')
loglog(modes,err_l2_deg_3_256,'-o','Linewidth',1.5,'DisplayName','p = 3, subs = 256')
loglog(modes,err_l2_deg_4_128,'-o','Linewidth',1.5,'DisplayName','p = 4, subs = 128')
loglog(modes,err_l2_deg_4_256,'-o','Linewidth',1.5,'DisplayName','p = 4, subs = 256')
legend('Location','southeast')
title('Error stability for solution of truncated Fourier expansions')
xlabel('modes / order of truncation')
ylabel('|\Re(u)-\Re(u_h)|_{L^2}')

figure ('Units', 'pixels', 'Position', [150 200 1000 350])
hold on; grid on;
loglog(modes,err_h1_deg_3_128,'-o','Linewidth',1.5,'DisplayName','p = 3, subs = 128')
loglog(modes,err_h1_deg_3_256,'-o','Linewidth',1.5,'DisplayName','p = 3, subs = 256')
loglog(modes,err_h1_deg_4_128,'-o','Linewidth',1.5,'DisplayName','p = 4, subs = 128')
loglog(modes,err_h1_deg_4_256,'-o','Linewidth',1.5,'DisplayName','p = 4, subs = 256')
legend('Location','southeast')
title('Error stability for solution of truncated Fourier expansions')
xlabel('modes / order of truncation')
ylabel('|\Re(u)-\Re(u_h)|_{H^1}')

