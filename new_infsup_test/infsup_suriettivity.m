function [mu, W, D] = infsup_suriettivity (A, Mv, Mw)
% infsup computes the discrete infsup constant for suriettivity, 
% i.e. the positive inf of 
% inf_w sup_v (A*v,w)/ norm(v) norm(w)
% by solving an eigenvalue problem
% A is the operator matrix, Mv and Mw are the
% bilinear form giving the norm

%[W, D] = eigs ( A*lsqminnorm(Mv,(A')), Mw, 1, 'smallestabs');
[W, D] = eigs ( A*(Mv\(A')), Mw, 1, 'smallestabs');
mu = sqrt(real(D));
end