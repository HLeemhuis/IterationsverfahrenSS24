% QCD Matrix test
load('conf5_0-4x4-26.mat');
kappa_c = 0.21070;

%load('conf6_0-8x8-30.mat');
%kappa_c = 0.15649;



D = Problem.A;
n = size(D,1);

kappa = 0.995*kappa_c;
A = speye(n)-kappa*D;

gamma_5 = [ 0 0 1 0; 0 0 0 1; 1 0 0 0; 0 1 0 0];
Gamma_5 = kron(speye(n/12), kron(gamma_5, speye(3)) );

Q = Gamma_5*A;
norm(Q-Q','fro')

%Q is Hermitian and indefinite, A has spectrum in the right half plane

rng(42);  % initialize random nunber generator
b = rand(n,1);

% stopping parameters
tol = 1e-10;
maxit  = 1000;
restart = 64;

%options.type = 'crout';
%options.milu = 'row';
%options.droptol = 0.1;
[L,U] = ilu(A);

%minres
%[x,flag,relres,iter,resvec_minres] = minres(Q,Gamma_5*b,tol,maxit);

%gmres
[x,flag,relres,iter,resvec_gmres] = gmres(A,b,restart,tol,maxit); 
[x,flag,relres,iter,resvec_pgmres] = gmres(A,b,restart,tol,maxit,L,U);
close all
figure
semilogy(resvec_gmres/norm(b));
hold on
semilogy(resvec_pgmres/norm(b));
legend('GMRES','prec. GMRES');
