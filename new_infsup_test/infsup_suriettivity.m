function [mu, W, D] = infsup_suriettivity (A, Mv, Mw)
% infsup computes the discrete infsup constant for suriettivity, 
% i.e. the positive inf of 
% inf_w sup_v (A*v,w)/ norm(v) norm(w)
% by solving an eigenvalue problem
% A is the operator matrix, Mv and Mw are the
% bilinear form giving the norm

%[W, D] = eigs ( A*(Mv\(A')), Mw, 1, 'smallestabs'); 

% see documentation for Afun in eigs... should use Afun = @(x) A\x !
% here we have (A * B * C )\x. This writes C \ (B \( A \x)).
% in our case C = A' and B = Mv^-1. Instead of computing the inverse and
% inverting again, we sould use the following: A' \ (Mv * ( A\ x) )
% therefore we use the following:
F = @(u) A' \ (Mv * (A\u)); 
[W, D] = eigs ( F, size(A,1), Mw, 1, 'smallestabs', 'IsFunctionSymmetric', true);
mu = sqrt(real(D));
end