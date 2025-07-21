function [q,beta] = reorth(v,W,Q)

for count = 1 : 2,
    v = v - [W,Q]*([W,Q]'*v);
end;

% % % % % % % % deflation % % % % % % % % 
tol = 1e-12;
[U,S,V] = svd(v,'econ');
i = length(diag(S));
while (i>0)&&(S(i,i)<tol)
    i = i-1;
end
if i==0
    i = length(diag(S));
end;
q = U(:,1:i);
beta = S(1:i,1:i)*V(:,1:i)';    
    
end

