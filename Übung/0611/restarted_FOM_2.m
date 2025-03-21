function [x] = restarted_FOM_2(x0, A, b, k, mr, tol)
% x0 initial guess, Matrix A, RHS b, restart-length k, maximum number of
% restarts mr, desired prec tol


r = b - A*x0;
beta = norm(r);
v = r /beta;

restarts = 0;
while restarts < mr && norm(r) > tol
    restarts = restarts + 1;
    iter = 0;
    while norm(r) > tol && iter < k
        iter = iter +1;
        e1 = zeros(iter, 1);
        e1(1) = 1;
    
        [V, H] = arnoldi(A, v, iter);
        y = beta *inv(H(1:iter, 1:iter)) * e1;
        x = x0 + V(:,1:iter) * y;
        r = b -  A * x;
    end
    v = V(:, iter+1);
    beta = - H(iter+1,iter) * y(k);
end