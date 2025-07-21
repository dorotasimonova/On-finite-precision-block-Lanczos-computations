function [T, V, beta, vlast, q] = lanczos(A, v, k, reo)

beta    = 0.0; 
n = size(A,1);
q = v/norm(v); 
v = zeros(n,1);

for j = 1 : k,
    
    V(:,j) = q;
    
    if j>1,
        T(j-1,j) = beta;
        T(j,j-1) = beta;
    end
        
    v = A*q - beta*v;
    alpha = q'*v;
    v = v - alpha*q;
    
    if j==k vlast=q; end;
    for count = 1 : reo,
        v = v - V*(V'*v);
    end;        
    beta = norm(v);
    q = v/beta;
    v = V(:,j);

    T(j,j) = alpha;    
    
end;

end
