function [sol, res] = tschebyscheff(A, b, x0, f1, f2, maxiter, tol)
%A als Matrix
%b RHS
%x0 Startvektor
%f1, f2 Fokuspunkte der Ellipse
%maxiter und toll sind selbsterklärend
%sol enthält Lösung
%res als array mit normen der residuen
af = 0.5*(f2 - f1);
bf = 0.5*(f2 + f1);

xlast = x0;
rlast = b - A*x0;
alphacurr = -2 / bf; 

xcurr = x0 - alphacurr * rlast;
rcurr = alphacurr * A * rlast  + rlast;

clast = 1;
ccurr = -bf /af;

res = zeros(1,maxiter);
iter = 1;
while (iter < maxiter && norm(rcurr) > tol)
    cnext = -2 /af * bf * ccurr - clast;
    alphanext = ccurr / cnext * 2 /af;
    xnext = - alphanext * bf * xcurr + (1 + alphanext * bf) * xlast  - alphanext * rcurr;
    rnext = b - A * xnext;

    clast = ccurr;
    ccurr = cnext;
    
    xlast = xcurr;
    xcurr = xnext;

    rcurr = rnext;

    res(iter) = norm(rcurr);
    iter = iter +1;
end

res = res(1:iter);
sol = xcurr;
