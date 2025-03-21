function [V,T] = lanczos(A,r,m)
% performs m steps of the Lanczos process to
% compute the Lanczos vectors which form an ONB 
% for the Krylov subspace K_m(A,r) 
% The Lanczos vectors are the columns of the nx(m+1)
% matrix V, and T is (m+1)xm, tridiagonal and contains
% the orthogonalization coefficients. The satisfy
% AV(:,1:m) = V(:,1:m+1)*T

end