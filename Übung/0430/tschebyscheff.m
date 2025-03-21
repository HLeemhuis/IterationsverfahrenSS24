function [sol, res] = tschebyscheff(x0, A, b, f1, f2, maxiter, tol)
% initial guess x0, Matrix A, RHS b, focuspoints of ellipse f1 and f2
% maximum number of iterations maxiter,
% maximum precision tol
% solution sol, array of residuum norm res 
af = 0.5 * (f2-f1);
bf = 0.5 * (f2+f1);

clast = 1;
ccurr = - bf / af;
alphacurr = - 1/bf;


xlast = x0;
rlast = b - A*x0;


xcurr = xlast - alphacurr * rlast;
rcurr = alphacurr * A * rlast + rlast;

iter = 1; %counter for iterations
res = zeros(1,maxiter); %contains norms of residuals
while(norm(rcurr) > tol && iter < maxiter)
    cnext = - ccurr  * 2 * bf/af - clast;
    alphanext = ccurr / cnext * 2 / af;

    xnext = - alphanext* bf * xcurr  + (1 + alphanext * bf) * xlast - alphanext * rcurr;
    rnext = b - A * xnext;


    xlast = xcurr;
    xcurr = xnext;
    clast = ccurr;
    rlast = rcurr;
    rcurr = rnext;
    ccurr = cnext;

    res(iter) = norm(rcurr); % store norm of residuum
    iter = iter+1;
end

res = res(1:iter-1);
sol = xnext;
end