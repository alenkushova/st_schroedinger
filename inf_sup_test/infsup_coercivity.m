function [mu, D] = infsup_coercivity(A,Mv,Mw)
% infsup computes the discrete infsup constant for coercivity, 
% i.e. the positive inf of 
% inf_v sup_w (A*v,w)/ norm(v) norm(w)
% by solving an eigenvalue problem
% A is the operator matrix, Mv and Mw are the
% bilinear form giving the norm
[~, D] = eigs ( A'*(Mw\A), Mv, 1, 'smallestabs');
mu = sqrt(real(D));
end