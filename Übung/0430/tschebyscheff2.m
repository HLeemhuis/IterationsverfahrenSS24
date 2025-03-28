function [sol, res] = tschebyscheff2(x0, A, b, f1, f2, maxiter, tol)
% initial guess x0, Matrix A, RHS b, focuspoints of ellipse f1 and f2
% maximum number of iterations maxiter,
% maximum precision tol
% solution sol, array of residuum norm res 
af = 0.5 * (f2-f1);
bf = 0.5 * (f2+f1);

dcurr = -bf / af;
xlast = x0;
rlast = b - A*xlast;

alphacurr = -2 / bf;
xcurr = xlast - alphacurr * rlast;
rcurr = alphacurr * A * rlast + rlast;

iter = 1; %counter for iterations
res = zeros(1,maxiter); %contains norms of residuals
while(norm(rcurr) > tol && iter < maxiter)
    dnext = -2 / af * bf - 1/dcurr;
    alphanext = 2 / (af*dnext);
    xnext = (-alphanext * bf) *xcurr + (1 + alphanext*bf) *xlast - alphanext*rcurr;
    rnext = b - A*xnext;
    

    xlast = xcurr;
    xcurr = xnext;
    
    rlast = rcurr;
    rcurr = rnext;

    dcurr = dnext;

    res(iter) = norm(rcurr); % store norm of residuum
    iter = iter+1;
end


res = res(1:iter-1);
sol = xnext;
end