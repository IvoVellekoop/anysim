function sims = default_simulations()
    % DEFAULT_SIMULATIONS() returns a struct array with
    % options for commonly used linear solvers.
    % This array can be used in compare_simulations
    sims = add([], "GMRES20", ...
        @(A, b, tol, Nit) gmres(A, b, 20, tol, Nit), itfactor = 20 + 1);

    sims = add(sims, "GMRES5", ...
        @(A, b, tol, Nit) gmres(A, b, 5, tol, Nit), itfactor = 5 + 1);

    sims = add(sims, "BiCGSTAB", ...
        @(A, b, tol, Nit) bicgstab(A, b, tol, Nit), itfactor = 2);

    sims = add(sims, "BiCGSTAB(l)",...
        @(A, b, tol, Nit) bicgstabl(A, b, tol, Nit), itfactor = 4);

    sims = add(sims, "Rich100", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 1.0));

    sims = add(sims, "Rich90", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.9));

    sims = add(sims, "Rich80", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.8));

    sims = add(sims, "Rich70", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.7));

    %sims = add(sims, "CGS",...
    %    @(A, b, tol, Nit) cgs(A, b, tol, Nit), itfactor = 2);
end

function simulations = add(simulations, name, fn, opts)
    arguments
        simulations
        name (1,1) string
        fn (1,1) function_handle
        opts.itfactor (1,1) double = 1
    end
    opts.name = name;
    opts.function = fn;
    simulations = [simulations opts];
end
