function sims = default_simulations(type, opt)
    arguments
        type (1,1) string {mustBeMember(type, ["symmetric", "nonsymmetric"])} = "nonsymmetric"
        opt.has_adjoint = false
        opt.N_restart_gmres (1,1) double = 5
    end

    % DEFAULT_SIMULATIONS(OPT) returns a struct array with
    % options for commonly used linear solvers.
    % This array can be used in compare_simulations
    sims = add([], "Rich. (α = 1.0)", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 1.0), ...
        symmetric = true, nonsymmetric = true);

    sims = add(sims, "Rich. (α = 0.75)", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.75), ...
        symmetric = true, nonsymmetric = true);

    sims = add(sims, "Rich. (α = 0.5)", ...
        @(A, b, tol, Nit) rich(A, b, tol, Nit, 0.5), ...
        symmetric = true, nonsymmetric = true);

    sims = add(sims, "GMRES (restart = " + opt.N_restart_gmres+")", ...
        @(A, b, tol, Nit) gmres(A, b, opt.N_restart_gmres, tol, Nit), ...
        symmetric = true, nonsymmetric = true, itfactor = opt.N_restart_gmres + 1);

    sims = add(sims, "GMRES (restart = " + opt.N_restart_gmres * 4+")", ...
        @(A, b, tol, Nit) gmres(A, b, opt.N_restart_gmres * 4, tol, Nit), ...
        symmetric = true, nonsymmetric = true, itfactor = opt.N_restart_gmres * 4 + 1);

    sims = add(sims, "BiCGSTAB",...
        @(A, b, tol, Nit) bicgstab(A, b, tol, Nit),...
        nonsymmetric = true, itfactor = 2);

    %sims = add(sims, "bicgstabl",...
    %    @(A, b, tol, Nit) bicgstabl(A, b, tol, Nit),...
    %    nonsymmetric = true, itfactor = 4);

    %sims = add(sims, "bicg",...
    %    @(A, b, tol, Nit) bicg(A, b, tol, Nit),...
    %    nonsymmetric = true, itfactor = 2, require_adjoint = true);

    %sims = add(sims, "qmr", ...
    %    @(A, b, tol, Nit) qmr(A, b, tol, Nit), ...
    %    nonsymmetric = true, itfactor = 2, require_adjoint = true);

%     sims = add(sims, "lsqr", ...
%         @(A, b, tol, Nit) lsqr(A, b, tol, Nit), ...
%         nonsymmetric = true, itfactor = 2, require_adjoint = true);
%     
    sims = add(sims, "CGS", ...
        @(A, b, tol, Nit) cgs(A, b, tol, Nit), ...
        nonsymmetric = true, symmetric = true, itfactor = 2);
% 
%     sims = add(sims, "pcg", ...
%         @(A, b, tol, Nit) pcg(A, b, tol, Nit), ...
%         symmetric = true, itfactor = 1);

    switch type
        case "symmetric"
            sims = sims([sims.symmetric]);
        case "nonsymmetric"
            sims = sims([sims.nonsymmetric]);
    end
    if ~opt.has_adjoint
        sims = sims(~[sims.require_adjoint]);
    end
end

function simulations = add(simulations, name, fn, opts)
    arguments
        simulations
        name (1,1) string
        fn (1,1) function_handle
        opts.itfactor (1,1) double = 1
        opts.require_adjoint (1,1) logical = false
        opts.symmetric (1,1) logical = false         % method should be tried for symmetric systems
        opts.nonsymmetric (1,1) logical = false      % method should be tried for nonsymmetric systems
    end
    opts.name = name;
    opts.function = fn;
    simulations = [simulations opts];
end
