function [V, T] = lanczos(A, r, m)

n = size(A,1);
V = zeros(n,m+1);
V(:,1) = r / norm(r);
T = zeros(m+1, m);

for i = 2:m+1
    V(:,i) = A * V(:,i-1);
    if i == 2
        V(:,i) = V(:,i);
    else
        V(:,i) = V(:,i) - T(i-1,i-2)*V(:,i-2); 
    end
    T(i-1,i-1) = V(:,i)' * V(:,i-1);
    V(:,i) = V(:,i) - T(i-1,i-1) * V(:,i-1);
    T(i,i-1) = norm(V(:,i));
    T(i-1,i) = T(i,i-1);
    V(:,i) = V(:,i) / T(i,i-1);
end
end