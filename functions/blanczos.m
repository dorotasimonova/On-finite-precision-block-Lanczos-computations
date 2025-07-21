function [T, V, vlast, beta, q] = blanczos(A, v, k)

n = size(A,1);
p = size(v,2);

beta = zeros(p); 
[q,~] = qr(v,0);
bj = 1:p;
v = zeros(n,p);

for j = 1 : k,
    
    V(:,bj) = q;
    
    if j>1,
        T(bj-p,bj) = beta';
        T(bj,bj-p) = beta;
    end
    
    v = A*q - v*beta';
    alpha = q'*v;
    v = v - q*alpha;
    
    if j==k vlast=q; end;        
    [q,beta] = qr(v,0);
    v = V(:,bj);
    
    T(bj,bj) = alpha;    
    bj = bj + p;
end;


