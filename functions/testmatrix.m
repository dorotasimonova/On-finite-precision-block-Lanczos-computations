function [A,b] = testmatrix(T,omega,p)

n = size(T,2)*p;
[U,~] = qr(randn(n,n));
A = U*kron(T,eye(p)+randn(p)*omega)*U';
n = size(T,2);
b = U*kron(eye(n,1),eye(p));

end

