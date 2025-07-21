function [Sm] = selectSm(T, beta, tol, p)

[S,theta]=eigsvd(T);
Sm = [];

for j = 1:length(theta)
    delta = max([eps,norm(beta*S(end-p+1:end,j))]);
    if (delta>tol)
        Sm = [Sm,S(:,j)];
    end
end

end

