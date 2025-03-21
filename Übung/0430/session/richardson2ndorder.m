function [sol, res] = richardson2ndorder(A, b, x0, omega, maxiter, tol)
%A als Matrix
%b RHS
%x0 Startvektor
%omega optimierungsparameter
%maxiter und toll sind selbsterklärend
%sol enthält Lösung
%res als array mit normen der residuen

xlast = x0;
xcurr = x0 - A * x0 + b;

res = zeros(1,maxiter);
iter = 1;
rcurr = b - A*xcurr;
while (iter < maxiter && norm(rcurr) > tol)
    xnext = omega * ( b + xcurr - A * xcurr) + ( 1 - omega)* xlast;

    rnext = b - A*xnext;

    xlast = xcurr;
    xcurr = xnext;
    rcurr = rnext;

    res(iter) = norm(rcurr);
    iter = iter +1;
end

res = res(1:iter);
sol = xcurr;
