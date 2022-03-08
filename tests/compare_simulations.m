function results = compare_simulations(sim, source, methods, opt)
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
%
%   OPT option structure with optional fields:
%       .analytical_solution    When present, also computes errors with
%               respect to analytical solution. Insert NaN for voxels to ignore.
%       .preconditioned

% note: simulation must have been done with opt.crop=false

defaults.analytical_solution = []; % no analytical solution given
defaults.tol = []; % auto: use residual of AnySim as tolerance.
defaults.iter = []; % auto: use same number of operator evaluations as AnySim
defaults.preconditioned = false;
opt = set_defaults(defaults, opt);
M = numel(methods);

%% First run the AnySim simulation
[u, state] = sim.exec(source);
results(M+1).name = 'AnySim';
results(M+1).value = gather(u);
results(M+1).iter = state.iteration;
fprintf('\nAnySim original: ');
sz = size(u);

if opt.preconditioned
    [A, state] = sim.preconditioned;
    b = sim.preconditioner(source);
    b = b(:);
    up = pagemtimes(inv(sim.Tr), u);     % u' = Tr^(-1) u
    results(M+1).residual = gather(norm(A(up(:))-b) / norm(b(:))); % Residue = ‖Γ(Ax-b)‖/‖b‖
else
    [A, state] = sim.operator;   % scaled operator L'+V'
    b = source(:);
    results(M+1).residual = gather(norm(A(u(:))-b(:)) / norm(b(:))); % Residue = ‖Ax-b‖/‖b‖
end


if isempty(opt.tol)
    tol = results(M+1).residual * gather(norm(b));
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
    state.iteration = 0;

    % run simulation and store results
    fprintf("\n" + m.name + ": ");
    [val, flag, relres, ~] = m.function(A, b, tol / norm(b), ceil(Nit / itfactor));
    results(m_i).flag = flag;
    results(m_i).iter = state.iteration;
    results(m_i).name = m.name;
    results(m_i).value = gather(pagemtimes(sim.Tr, reshape(val, sz))); % compensate for scaling of operator A
    
    % Relative residual =: ‖Ax-b‖/‖b‖
    results(m_i).residual = gather(relres); %norm(A(val)-b)/norm(b);
end

%% Compare with analytical theory
if ~isempty(opt.analytical_solution)
    mask = ~isnan(opt.analytical_solution(:));
    a = opt.analytical_solution(mask);
    for r_i = 1:length(results)
        err = results(r_i).value(mask) - a;
        results(r_i).rel_error = norm(err)/norm(a);
    end
end

result_summary = rmfield(results, 'value');
result_summary = rmfield(result_summary, 'flag');
disp(' ');
disp(struct2table(result_summary));


