function [V,T] = lanczos(A,r,m)
% performs m steps of the Lanczos process to
% compute the Lanczos vectors which form an ONB 
% for the Krylov subspace K_m(A,r) 
% The Lanczos vectors are the columns of the nx(m+1)
% matrix V, and T is (m+1)xm, tridiagonal and contains
% the orthogonalization coefficients. The satisfy
% AV(:,1:m) = V(:,1:m+1)*T

% Algorithm  6.6 

n = size(A,1);
V = zeros(n,m+1);
T = zeros(m+1,m);
V(:,1) = r/norm(r);
beta_old = 0;

for i = 2:m+1
    if i == 2
        q = A * V(:,i-1);
    else    
        q = A * V(:,i-1) - beta_old * V(:,i-2);
    end
    alpha_old = q' * V(:,i-1);
    q = q - alpha_old * V(:,i-1);

    beta_curr = norm(q);
    V(:,i) = q/beta_curr;
    beta_old = beta_curr;
    if i ~= m+1
        T(i-1,i-1) = alpha_old;
        T(i,i-1) = beta_old;
        T(i-1,i) = beta_old;
    end
end

T(m,m) = alpha_old;
T(m+1,m) = beta_old;

end