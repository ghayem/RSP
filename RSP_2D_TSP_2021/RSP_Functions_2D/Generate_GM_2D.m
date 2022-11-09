

function [GM] = Generate_GM(roh, sig, x, m, n)


C = My_2D_Cov ( roh, sig, x, x);

K = size(x,1);

% [s v d] = svd(C);
% EPS = max(v(:))/2;
% L = chol(C + (EPS)*eye(K), 'lower');

L = chol(C + (1e-6)*eye(size(C)), 'lower');

% L = chol(C , 'lower');
m1 = repmat(m,1,n);
GM = L * randn(K,n) + m1;

end