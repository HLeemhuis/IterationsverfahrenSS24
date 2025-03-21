function [V, H] = arnoldi_RE(A,r,m)
n = size(A,1);
V = zeros(n,m+1);
H = zeros(m+1,m);
V(:,1) = r/norm(r);

for i = 2:m+1
    V(:,i) = A * V(:,i-1);
    for j = 1:i-1
        H(j,i-1) = V(:,i)' * V(:,j);
        V(:,i) = V(:,i) - H(j,i-1) * V(:,j);
    end
    H(i,i-1) = norm(V(:,i));
    V(:,i) = V(:,i) / H(i,i-1);
end

for i = 2:m+1
    for j = 1 : i-1
        alpha = V(:,i)'*V(:,j);
        H(j,i-1) = alpha + H(j,i-1);
        V(:,i) = V(:,i) - alpha * V(:,j);
    end
    V(:,i) = V(:,i) / norm(V(:,i));
end