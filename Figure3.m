
addpath('functions')

clear
close all 

% % % % % % % % initialization % % % % % % % %
p=2;
n=48;   
rho = 0.8;
lambda1 = .001; lambdan = 1;        % % % % ||A|| = 1
% lambda1 = .1; lambdan = 100;        % % % % ||A|| = 100
[Asv, ~] = strakosmatrix(n, rho, lambda1, lambdan);
[Tsv, ~] = lanczos(Asv, randn(n,1), n, 2);
omega = 1e-12;
[A,b] = testmatrix(Tsv, omega, p);
n = size(A,2);

% % % % % % % % computations % % % % % % % % 
k = 50;
[T, V, ~, beta, v] = blanczos(A, b, k);
tol = 1e-5*norm(A); 
V = [V,v];
bi = [1:p];
T(end+bi,flip(end-bi+1)) = beta; 
for i = 1:k-1
    bi = bi+p;
    [Sm] = selectSm(T(1:bi(end),1:bi(end)), T(bi+p,bi), tol, p);
    if isempty(Sm)
        fprintf('No Ritz vector fulfils the criterion for Zm.\n')
        return
    end;
    Zm = V(:,1:bi(end))*Sm; [W,R] = qr(Zm,0);
    rk = [V(:,bi)'*V(:,1:(bi(end)-p)),zeros(p)];
    norms1(i) = norm(W'*V(:,bi+p)*T(bi+p,bi));
    norms2(i) = norm(T(bi+p,bi)*rk*Sm*inv(R)); 
    norms3(i) = norm((eye(size(W,1)) - W*W')*(V(:,bi)*T(bi+p,bi)'));
end

% % % % % % % % plots % % % % % % % % 
figure(1)
semilogy(norms1,'b*-'); hold on
semilogy(norms2,'ro-'); hold on
semilogy(norms3,'kx-'); hold on
yline(tol,'--');
legend('||W_k^Tv_{k+1}\beta_{k+1}||','||\beta_{k+1}r_k^TS_mR_k^{-1}||','||(I - W_kW_k^T)v_k\beta^T_{k+1}||')
