function [results, table] = compare_simulations(sim, source, methods, opt)
    arguments
        sim AnySim
        source
        methods
        opt.analytical_solution = []; % no analytical solution given
        opt.tol = 1E-3; % relative residual ‖Ax-b‖/‖b‖ termination condition. [] = auto: use residual of AnySim as tolerance.
        opt.iter = 25000; % [] = auto: use same number of operator evaluations as AnySim
        opt.preconditioned logical = true;
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
    
    if opt.preconditioned
        [A, state] = sim.preconditioned;
        b = sim.preconditioner(source);
        state.source = b;
        b = b(:);
    else
        [A, state] = sim.operator;   % scaled operator L'+V'
        b = source(:);
        state.source = b;
        up = pagemtimes(inv(sim.Tr), u);     % u' = Tr^(-1) u
        results(M+1).residual = gather(norm(A(up(:))-b(:)) / norm(b(:))); % Residue = ‖Ax-b‖/‖b‖
    end


if isempty(opt.tol)
    tol = results(M+1).residual;
else
    tol = opt.tol;
end

% Number of iterations (operator evaluations actually)
if isempty(opt.iter)
    Nit = results(M+1).iter;
else
    Nit = opt.iter;
end

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
    % reset interation count
    state.reset();

    % run simulation and store results
    fprintf("\n" + m.name + ": ");
    [val, flag, relres, ~] = m.function(A, b, tol, max(ceil(Nit / itfactor - 1), 1));
    state.finalize();
    results(m_i).flag = flag;
    results(m_i).iter = state.iteration + 1; %1 extra for computing preconditioned sources
    results(m_i).name = m.name;
    results(m_i).value = sim.finalize(val);
    results(m_i).time = state.run_time;
        
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

rstore = results(1:end-1);
disp(' ');
disp(struct2table(result_summary));
S = dbstack();
filename = S(2).file(6:end-2); % remove test_ and .m
disp(filename)

header = "";
for r = rstore
    name = strrep(strrep(r.name, 'α', '$\alpha$'), '_', ' ');
    header = header + sprintf("& \\rotatebox{90}{%s}", name);
end
header = header + "\\";

data = sprintf("%s ", strrep(filename, '_', ' '));
for r = rstore
    if r.iter < Nit && r.flag == 0
        data = data + sprintf("& %d", r.iter);
    elseif r.flag == 3
        data = data + sprintf("& s"); % stagnates, may be due too low machine precision
    else
        data = data + sprintf("& -");
    end
end
data = data + "\\";
table = [header, data];
f = fopen([filename '.log'], 'w+');
fprintf(f, "%s\n", header);
fprintf(f, "%s\n", data);
fclose(f);
