function results = compare_simulations(sim, source, methods, opt)
    arguments
        sim AnySim
        source
        methods
        opt.analytical_solution = []; % no analytical solution given
        opt.measure_time logical = false;
        opt.nb_time_measurements = 1;
    end
    %% Helper function to compare different simulation algorithms
    %   COMPARE_SIMULATIONS(SIM, SOURCE, METHODS, OPT) 
    %   Executes each of the iterative methods in METHODS for the given
    %   simulation object SIM and given SOURCE.
    %
    %   METHODS is a struct array, with each entry containing a .name
    %   with the display name of the method, and a .function containing a
    %   function handle of the type @(A, b, Nit), taking three parameters
    %   (operator A, source b, number of iterations Nit).
    %   For example:
    %       simulations(1).name = 'gmres';
    %       simulations(1).function = @(A, b, Nit) gmres(A, b, Nrestart_gmres, 1E-10, Nit);
    %       compare_simulations(sim, source, simulations)

    M = numel(methods);
    
    %% First run the AnySim simulation
    fprintf('\nAnySim original: ');
    [u, state] = sim.exec(source);
    results(M+1).name = 'AnySim';
    results(M+1).value = gather(u);
    results(M+1).iter = state.iteration;
    results(M+1).residual = state.residuals(end);
    results(M+1).time = state.run_time;
    results(M+1).flag = 0;
    
    [A, state] = sim.preconditioned;
    b = sim.preconditioner(source, state);
    state.source = b;
    b = b(:);
    
    tol = sim.opt.termination_condition.tolerance;
    Nit = sim.opt.termination_condition.iteration_count;
    
    %% Run all other simulations
    for m_i = 1:M
        m = methods(m_i);
        % Determine the maximum number of iterations
        % compensate for the fact that some methods
        % perform multiple evaluations of 'A' per iteration
        if isfield(m, 'itfactor')
            itfactor = m.itfactor;
        else
            itfactor = 1;
        end
        % Repeat each calculation multiple times and make a table to compare the minimum execution times
        if opt.measure_time
            nb_time_measurements = opt.nb_time_measurements;
        else
            nb_time_measurements = 1;
        end
        run_times = zeros(1, nb_time_measurements);
        for measurement_idx = 1:nb_time_measurements
            % reset interation count
            state.reset();
    
            % run simulation and store results
            fprintf("\n");
            if opt.measure_time
                fprintf("\n time measurement %d/%d ", measurement_idx, nb_time_measurements);
            end
            fprintf(m.name + ": ");
            try
                [val, flag, relres, ~] = m.function(A, b, tol, max(ceil(Nit / itfactor - 1), 1));
            catch err
                warning(err.message);
            end
            state.finalize();
            run_times(measurement_idx) = state.run_time;
        end
        results(m_i).flag = flag;
        results(m_i).iter = state.iteration + 1; %1 extra for computing preconditioned sources
        results(m_i).name = m.name;
        results(m_i).value = sim.finalize(val);
        results(m_i).time = min(run_times); % state.run_time;
            
        % Relative residual =: ‖Ax-b‖/‖b‖
        results(m_i).residual = gather(relres); %norm(A(val)-b)/norm(b);
    end
    
    %% Compare with analytical theory
    if ~isempty(opt.analytical_solution)
        if size(opt.analytical_solution) ~= size(results(end).value)
            error("Size of the analytical solution does not match size of the simulation");
        end
        mask = ~isnan(opt.analytical_solution(:));
        a = opt.analytical_solution(mask);
        for r_i = 1:length(results)
            err = results(r_i).value(mask) - a;
            results(r_i).rel_error = norm(err)/norm(a);
        end
    end
    
    result_summary = rmfield(results, 'value');
    result_summary = rmfield(result_summary, 'flag');
    disp('\n');
    disp(struct2table(result_summary));
    