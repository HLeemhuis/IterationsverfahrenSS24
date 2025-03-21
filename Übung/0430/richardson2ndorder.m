function [sol, res] = richardson2ndorder(x0, A, b, alpha, maxiter, tol)
% initial guess x0, Matrix A, RHS b
% maximum number of iterations maxiter,
% maximum precision tol
% solution sol, array of residuum norm res 


xlast = x0;
xcurr = x0 - A*x0 + b;

iter = 1; %counter for iterations
res = zeros(1,maxiter); %contains norms of residuals
rnext = b - A*xcurr;

while (iter < maxiter && norm(rnext) > tol)
    xnext = alpha*(xcurr - A* xcurr + b) + (1- alpha) * xlast;
    rnext = b - A*xnext;

    xlast = xcurr;
    xcurr = xnext;

    res(iter) = norm(rnext); % store norm of residuum
    iter = iter+1;
end

res = res(1:iter-1);
sol = xnext;
end
