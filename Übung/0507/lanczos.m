function [V,T] = lanczos(A,r,m)

n = size(A,1);

V = zeros(n,m+2);
V(:,2) = r/norm(r);
beta_last = 0;
T = zeros(m+1,m);

for i = 2:m+1
    q = A * V(:,i) - beta_last*V(:,i-1);
    alpha_last = q' * V(:,i);
    q = q - alpha_last * V(:,i);
    beta_curr = norm(q);
    V(:,i+1) = q / beta_curr;
    T(i-1,i-1) = alpha_last;
    T(i, i-1) = beta_curr;
    if i <= m
        T(i-1,i) = beta_curr;
    end
    beta_last = beta_curr;
end
V = V(:,2:m+2);
end