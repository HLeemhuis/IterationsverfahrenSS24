function [x, relres] = RFOM_Session(A, b, x0, k, mr,tol)
%restart-länge k, maximale anzahl an restarts mr

r = b - A*x0;

beta = norm(r);

v = r/beta;

relres = zeros(1, k*mr);


restarts = 0;

while (norm(r) > tol) && restarts < mr
    restarts = restarts + 1;
    iter = 0;
    while (norm(r) > tol) && iter < k
        iter = iter + 1;
    
        e1 = zeros(iter, 1);
        e1(1) = 1;
    
        [V, H] = arnoldi(A, v, iter);
        vm1 = V(:,iter+1);
        hm1m = H(iter+1,iter);
        V = V(:,1:iter);
        H = H(1:iter,:);
    
        y = beta * inv(H) *e1;

        x = x0 + V*y;
    
        r = b- A*x;
        relres(iter + k*(restarts-1)) = norm(r);
    end
    v = vm1;
    beta = - hm1m * y(k);
end

relres = relres(1:iter + k*(restarts-1));