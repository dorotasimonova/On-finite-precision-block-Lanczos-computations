
addpath('functions')

clear
close all 

% % % % % % % % initialization % % % % % % % %
p=2;
n=48;   
rho = 0.8;
% lambda1 = .001; lambdan = 1;        % % % % ||A|| = 1
lambda1 = .1; lambdan = 100;        % % % % ||A|| = 100
[Asv, ~] = strakosmatrix(n, rho, lambda1, lambdan);
[Tsv, ~] = lanczos(Asv, randn(n,1), n, 2);
omega = 1e-12;
[A,b] = testmatrix(Tsv, omega, p);
n = size(A,2);

% % % % % % % % constructing T_N % % % % % % % %
k = 24; 
[T, V, ~, beta, ~] = blanczos(A, b, k);

tol = 1e-5*norm(A); % sqrt(k*n*p*eps)*norm(A); % 
[Sm] = selectSm(T, beta, tol, p);
if isempty(Sm)
    fprintf('No Ritz vector fulfils the criterion for Zm.')
    return
end;
Zm = V*Sm;
[W,R] = qr(Zm,0);

[T_N, normsHk] = contprocess(A, T, V(:,end-2*p+1:end-p), V(:,end-p+1:end), W);

% % % % % % % % computations % % % % % % % %
[~,lambda] = eigsvd(A);
[~,theta] = eigsvd(T_N);
clusters = zeros(2,length(lambda));
for i = 1:length(theta)
    [dist,ind] = min(abs(theta(i)-lambda));
    clusters(1,ind) = clusters(1,ind) + 1;
    if dist > clusters(2,ind) 
        clusters(2,ind) = dist; 
    end;
end;

% % % % % % % % plots % % % % % % % % 
figure(1)
semilogy([0:size(normsHk,2)-1], normsHk,'*-')
yline(tol,'--')

figure(2)
subplot(2,1,1);bar(log10(clusters(2,:)/(sqrt(eps)*norm(A))));
xlabel('\lambda_i'); ylabel('size of cluster / {\surd\epsilon} ||A||');
subplot(2,1,2);bar(clusters(1,:)); hold on;
xlabel('\lambda_i'); ylabel('# of Ritz vals in cluster');

