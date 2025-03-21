function A = konv_diff(N,eps)
% builds the matrix for the 5-point discretization
% of a convection-diffusion operator - Delta u + eps*u_x
% discretized on an N x N uniform grid in [0,1] x [0,1]

h = 1/(N+1);
e = ones(N,1);
C = spdiags([-e,e],[-1 1],N,N);
A = gallery('poisson',N) * h*eps*kron(speye(N),C);
end