function [x] = FOM(x0, A, b, tol, max_iter)
% x0 initial guess, Matrix A, RHS b, desired prec tol
r = b - A*x0;
beta = norm(r);
v = r /beta;

iter = 0;
while norm(r) > tol && iter < max_iter
    iter = iter +1;
    e1 = zeros(iter, 1);
    e1(1) = 1;

    [V, H] = arnoldi(A, v, iter);
    x = x0 + V(:,1:iter) * beta *inv(H(1:iter, 1:iter)) * e1;
    r = b -  A * x;
end