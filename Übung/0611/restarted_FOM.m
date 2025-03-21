function [x] = restarted_FOM(x0, A, b, k, mr, tol)
% x0 initial guess, Matrix A, RHS b, restart-length k, maximum number of
% restarts mr, desired prec tol



r = b - A*x0;
beta = norm(r);
v = r /beta;

iter = 0;
while norm(r) > tol && iter < mr
    iter = iter +1;

    e1 = zeros(k, 1);
    e1(1) = 1;

    [V, H] = arnoldi(A, v, k);
    hm1m = H(k+1,k);
    vm1 = V(:,k+1);
    V = V(:,1:k);
    H = H(1:k,1:k);
    
    y =  beta *inv(H) * e1;
    %r_k = beta * v_k+1
    beta = - hm1m * y(k);

    x = x0 + V * y;
    v = vm1;
    r = b -  A * x;
end