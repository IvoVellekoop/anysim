function simulations = default_simulations(varargin)
% DEFAULT_SIMULATIONS(OPT) returns a struct array with
% options for the most commonly used linear solvers.
% This array can be used in compare_simulations
    defaults.N_restart_gmres = 5;
    opt = set_defaults(defaults, varargin{:});
    simulations(1).name = 'gmres';
    simulations(1).itfactor = opt.N_restart_gmres + 1;
    simulations(1).function = @(A, b, tol, Nit) gmres(A, b, opt.N_restart_gmres, tol, Nit);
    simulations(2).name = 'bicgstab';
    simulations(2).itfactor = 2; % bicgstab performs two evaluations per iteration
    simulations(2).function = @(A, b, tol, Nit) bicgstab(A, b, tol, Nit);
    simulations(3).name = 'bicgstabl';
    simulations(3).itfactor = 4;
    simulations(3).function = @(A, b, tol, Nit) bicgstabl(A, b, tol, Nit);
    simulations(4).name = 'cgs';
    simulations(4).itfactor = 2;
    simulations(4).function = @(A, b, tol, Nit) cgs(A, b, tol, Nit);
end