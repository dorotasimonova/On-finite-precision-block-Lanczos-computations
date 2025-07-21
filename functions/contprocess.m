function [T, normsHk] = contprocess(A, T, v, q, W)

m = size(W,2); 
p = size(q,2);
n = size(A,2);
it = n-m;
bj = [size(T,2)-p+1:size(T,2)]; aj = bj;

beta = T(bj-p,bj)'; alpha = T(bj,bj);
v = A*q - q*alpha - v*beta';
[qnew,beta] = reorth(v,W,[]);
h = norm(v - qnew*beta);
v = q; q = qnew;

Q = []; normsHk = []; j = 1;

while (norm(double(beta)) > eps*norm(A)) && (j <= it) 
  
    Q = [Q,q];
    normsHk = [normsHk,h];
    
    aj = aj(end)+[1:size(beta,1)]; 
    T(bj,aj) = beta';
    T(aj,bj) = beta; 
    
    v = A*q - v*beta';
    alpha = q'*v; 
    v = v - q*alpha;
    [qnew,beta] = reorth(v,W,Q);
    h = norm(v - qnew*beta);
    
    v = q;
    q = qnew;
    
    T(aj,aj) = alpha;
    bj = aj; 
    j = j+1;

end

normsHk = [normsHk,h];

end

