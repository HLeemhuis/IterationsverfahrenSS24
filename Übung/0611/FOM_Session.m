function [x, relres] = FOM_Session(A, b, x0, tol, maxiter)

r = b - A*x0;

beta = norm(r);

v1 = r /beta;

relres = zeros(1, maxiter);


iter = 0;
while (norm(r) > tol) && iter < maxiter
    iter = iter + 1;

    e1 = zeros(iter, 1);
    e1(1) = 1;

    [V, H] = arnoldi(A, v1, iter);
    V = V(:,1:iter);
    H = H(1:iter,:);

    y = beta * inv(H) *e1;
    x = x0 + V*y;

    r = b- A*x;
    relres(iter) = norm(r);
end

relres = relres(1:iter);