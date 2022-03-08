function simulations = default_simulations(opt)
    arguments
        opt.N_restart_gmres (1,1) double = 5
    %    opt.alpha_rich (1,1) double = 0.75
    end
    % DEFAULT_SIMULATIONS(OPT) returns a struct array with
    % options for commonly used linear solvers.
    % This array can be used in compare_simulations
    simulations(1).name = "rich (α = 1.0)";
    simulations(1).function = @(A, b, tol, Nit) rich(A, b, tol, Nit, 1.0);
    simulations(1).itfactor = 1;
    simulations(2).name = "rich (α = 0.75)";
    simulations(2).function = @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.75);
    simulations(2).itfactor = 1;
    simulations(3).name = "rich (α = 0.5)";
    simulations(3).function = @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.5);
    simulations(3).itfactor = 1;
    simulations(6).name = "gmres";
    simulations(6).itfactor = opt.N_restart_gmres + 1;
    simulations(6).function = @(A, b, tol, Nit) gmres(A, b, opt.N_restart_gmres, tol, Nit);
    simulations(7).name = "bicgstab";
    simulations(7).itfactor = 2; % bicgstab performs two evaluations per iteration
    simulations(7).function = @(A, b, tol, Nit) bicgstab(A, b, tol, Nit);
    simulations(8).name = "bicgstabl";
    simulations(8).itfactor = 4;
    simulations(8).function = @(A, b, tol, Nit) bicgstabl(A, b, tol, Nit);
    simulations(4).name = "cgs";
    simulations(4).itfactor = 2;
    simulations(4).function = @(A, b, tol, Nit) cgs(A, b, tol, Nit);
    simulations(5).name = "pcg";
    simulations(5).itfactor = 2;
    simulations(5).function = @(A, b, tol, Nit) cgs(A, b, tol, Nit);
end