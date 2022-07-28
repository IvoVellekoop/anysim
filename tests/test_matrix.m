%% Experimental: tests solving a large matrix

kappa = 1000;
N = 1E5; 
density = 1E-3;
if ~exist('A', 'var') || ~exist('b', 'var')
    %A = gallery('randsvd', N, kappa, [], [], [], 1); % random matrix with condition number kappa
    %Nx = 10;
    %Ny = 10;
    %A = gallery('wathen', Nx, Ny);
    %N = size(A,1);
    %N = 1E3;
    %rows = repmat(1:N, 1, 10);
    %M = length(rows); % number of non-zero elements
    %A = sparse(rows, randi(N, M, 1), randn(M, 1), N, N);
    if exist("A.mat", "file")
        precomputed = load("A.mat"); % sprand(N, N, density, 1/kappa);
        A = precomputed.A;
        b = precomputed.b;
    else
        A = sprand(N, N, density, 1/kappa);
        b = randn(N, 1);
        save("A.mat", "A", "b", "density", "N", "kappa");
    end
end
xcorrect = lsqr(A, b, 1E-10, 200);

tol = 1E-8;
maxit = 1E3;
opt = AnySimOptions;
opt.precision = "double";
opt.termination_condition.relative_limit = tol;
opt.termination_condition.iteration_count = maxit;
opt.gpu_enabled = false;

opt.alpha = 1;
sim = MatrixSolve(A, opt);

simulations = default_simulations();
simulations(1:6) = []; % only execute AnySim original and conjugate gradient
[precond, table] = compare_simulations(sim, b, simulations, tol=tol, analytical_solution=xcorrect);
disp(norm(A * precond(1).value - b)/norm(b))
