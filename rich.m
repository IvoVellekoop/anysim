function [x,flag,relres,iter,resvec] = pcg(A,b,tol,maxit,alpha)
%RICH   Modified Richardson iteration, based on pcg.m
%   X = RICH(A,B) accepts a function handle A.
%   A(X) accepts a vector input X and returns the matrix-vector product
%   A*X.
%
%   X = RICH(A,B,TOL) specifies the tolerance of the method. If TOL is []
%   then RICH uses the default, 1e-6.
%
%   X = RICH(A,B,TOL,MAXIT) specifies the maximum number of iterations. If
%   MAXIT is [] then RICH uses the default, min(N,20).
%

%   [X,FLAG] = RICH(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    4 one of the scalar quantities calculated during RICH became too
%      small or too large to continue computing.
%
%   [X,FLAG,RELRES] = RICH(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If FLAG is 0, then RELRES <= TOL.
%
%   [X,FLAG,RELRES,ITER] = RICH(A,B,...) also returns the iteration number
%   at which X was computed: 0 <= ITER <= MAXIT.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = RICH(A,B,...) also returns a vector of the
%   estimated residual norms at each iteration including NORM(B-A*X0).
%
%   See also PCG, BICG, BICGSTAB, BICGSTABL, CGS, GMRES, LSQR, MINRES, QMR,
%   SYMMLQ, TFQMR, ICHOL, FUNCTION_HANDLE.

%   Copyright 1984-2020 The MathWorks, Inc.
%   Copyright 2022 Ivo Vellekoop
% 
if (nargin < 2)
    error(message('MATLAB:pcg:NotEnoughInputs'));
end

if ~iscolumn(b)
    error(message('MATLAB:pcg:RSHnotColumn'));
end
n = numel(b);

% Assign default values to unspecified parameters
if (nargin < 3) || isempty(tol)
    tol = 1e-6;
end

if tol <= eps
    warning(message('MATLAB:pcg:tooSmallTolerance'));
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:pcg:tooBigTolerance'));
    tol = 1-eps;
end
if (nargin < 4) || isempty(maxit) || maxit < 0
    maxit = min(n,20);
end
if (nargin < 5)
    alpha = 0.75;
end
if (nargin > 5)
    error(message('MATLAB:pcg:TooManyInputs'));
end

% Set up for the method
x = zeros(n,1); 
n2b = norm(b);
flag = 1;

resvec = zeros(maxit+1,1);         % Preallocate vector for norm of residuals
resvec(1,:) = n2b;                 % resvec(1) = norm(b-A*x0)

% loop over maxit iterations (unless convergence or failure)
for iter = 1 : maxit
    r = A(x) - b;   % residual
    x = x - alpha * r;
    normr = norm(r);
    relres = normr / n2b;
    resvec(iter+1,1) = normr;
        
    if ((normr == 0) || isinf(normr))
        flag = 4;
        break
    end

    % check for convergence
    if (relres <= tol)
        flag = 0;
        break
    end
end

resvec = resvec(1:iter+1,:);

