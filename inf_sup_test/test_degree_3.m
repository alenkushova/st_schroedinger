p = 3; 
eigenv = zeros(1,5);
infsup = zeros(1,5);
h = zeros(1,5);
for i = 2 : 6 
   nel = 2^i;
   [mu, ei]  = infsuptest (p, nel);
   infsup (i-1) = mu;
   eigenv (i-1) = ei;
   h (i-1) = 1/nel;
end
filename = ['degree_' num2str(p) '_infsup_results.mat'];
save( filename , 'h', 'infsup', 'eigenv')
fprintf (['Results saved in file: ' filename '\n\n'])