function [A, lambda, X] = strakosmatrix(n, rho, lambda1, lambdan)

lambda = lambda1*ones(n,1);
for i=2:n, lambda(i) = lambda(1) + ((i-1)/(n-1))*(lambdan-lambda1)*rho^(n-i); end;
% lambda(n-1) = lambda(n);
[X,~] = qr(randn(n,n));
A = X*diag(lambda)*X';

end

