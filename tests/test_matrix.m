kappa = 10;
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
        xcorrect = precomputed.xcorrect;
    else
        A = sprand(N, N, density, 1/kappa);
        b = randn(N, 1);
        xcorrect = bicgstab(A, b, 1E-8);
        save("A.mat", "A", "b", "xcorrect", "density", "N", "kappa");
    end
end
tol = 1E-8;
opt = AnySimOptions;
opt.precision = "double";
opt.termination_condition.relative_limit = tol;
opt.termination_condition.iteration_count = 1E3;
opt.alpha = 1;
sim = MatrixSolve(A, opt);

simulations = default_simulations("symmetric");
simulations(1:6) = []; % only execute AnySim original and conjugate gradient
precond = compare_simulations(sim, b, simulations, tol=tol);
