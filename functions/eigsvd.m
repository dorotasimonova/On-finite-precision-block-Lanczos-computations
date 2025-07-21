function [U,s] = eigsvd(A)
      
[U,S,V] = svd(A);
U = flip(U,2); V = flip(V,2);
U = 1/2*(U + V);
s = flip(diag(S));

end

