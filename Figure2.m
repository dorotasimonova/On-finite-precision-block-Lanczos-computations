
addpath('functions')

clear
close all 

% % % % % % % % initialization % % % % % % % %
p=2;
n=48;   
rho = 0.8;
lambda1 = .001; lambdan = 1;        % % % % ||A|| = 1
% lambda1 = .1; lambdan = 100;        % % % % ||A|| = 100
[A, lambda] = strakosmatrix(n, rho, lambda1, lambdan);
b = randn(n,p);

% % % % % % % % computations % % % % % % % %  
k = 60;
[T, V, ~, beta, v] = blanczos(A, b, k);
bi = [1:p];
V = [V,v];
T(end+bi,flip(end-bi+1)) = beta; 
for i = 1:k
    Fk = A*V(:,bi) - V*T(:,bi);
    norms1(i) = norm(Fk);
    norms2(i) = norm(V(:,bi)'*V(:,bi+p)*T(bi+p,bi));
    norms3(i) = norm(V(:,bi)'*V(:,bi)-eye(p));
    bi = bi + p;
end;

% % % % % % % % plots % % % % % % % % 
figure(1)
semilogy(norms1,'b-*'); hold on
semilogy(norms2,'r-x'); hold on
semilogy(norms3,'k-o');
yline(n*p*eps,'--');
yline(n*p*eps*norm(A),':');
legend('||\Delta v_j||','||v_j^Tv_{j+1}\beta_{j+1}||','||v_j^Tv_j-I||')
