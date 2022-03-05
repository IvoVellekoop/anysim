function results = compare_simulations(sim, source, methods, varargin)
%% Helper function to compare different simulation algorithms
%   COMPARE_SIMULATIONS(SIM, SOURCE, SIMULATIONS) executes AnySim
%   simulation SIM with the given SOURCE. Then executes each of the
%   competing methods in SIMULATIONS (see below). 
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

defaults.analytical_solution = [];
defaults.tol = []; % auto: use residual of AnySim as tolerance.
defaults.iter = []; % auto: use same number of operator evaluations as AnySim
defaults.preconditioned = false;
opt = set_defaults(defaults, varargin{:});



%% First run the AnySim simulation
[u, state] = sim.exec(source);

results(1).name = 'AnySim';
results(1).value = gather(u);
results(1).iter = state.iteration;
sz = size(u);

if opt.preconditioned
    [A, state] = sim.preconditioned;   % scaled operator L'+V'
    b = sim.preconditioner(source);      % returns s' = Tl s
else
    [A, state] = sim.operator;   % scaled operator L'+V'
    b = source; %.to_array();      % returns s' = Tl s
end
b = b(:);

% Residue = ‖Ax-b‖
up = pagemtimes(inv(sim.Tr), u); % u' = Tr^(-1) u
results(1).residual = norm(A(up(:))-b) / norm(b);

if isempty(opt.tol)
    tol = results(1).residual * norm(b);
else
    tol = opt.tol;
end

% Number of iterations (operator evaluations actually)
if isempty(opt.iter)
    Nit = state.iteration;
else
    Nit = opt.iter;
end

%% Run all other simulations
for m_i = 1:length(methods)
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
    [val, flag, relres, iter] = m.function(A, b, tol / norm(b), ceil(Nit / itfactor));
    results(m_i+1).flag = flag;
    results(m_i+1).iter = state.iteration;
    results(m_i+1).name = m.name;
    results(m_i+1).value = gather(pagemtimes(sim.Tr, reshape(val, sz))); % compensate for scaling of operator A
    
    % Relative residual =: ‖Ax-b‖/‖b‖
    results(m_i+1).residual = norm(A(val)-b)/norm(b);
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
disp(struct2table(result_summary));
% 
% plot(z, Iz, 'LineWidth', 1.5);
% hold on;
% 
% 
% ylimits = [0, 10];
% plot(z_inside, I_th, '--', 'LineWidth', 2.5);
% plot([zl, zl], ylimits, 'k:', 'LineWidth', 1); 
% plot([zr, zr], ylimits, 'k:', 'LineWidth', 1); 
% ylim(ylimits);
% xlabel('z [\mum]');
% ylabel('u [J/m^3]');
% set(gca, 'FontSize', 14);
% legend('simulated', 'theoretical', 'sample boundaries', 'FontSize', 12);
% title('Intensity profile inside diffusing slab')
% %
% relerr = norm(Iz(z_indices).'-I_th(:))/norm(I_th(:));
% disp(['Relative error (analytical) ' num2str(relerr)]);
% 
% %% Determine residual ||A x - b|| / ||b||
% 
% relerrs = norm(A(up(:))-b) / norm(b);
% disp(['Relative residual (source) ' num2str(relerrs)]);
% 
% %% Solve the problem using standard algorithms
% 
% %% GMRES
% disp('Running gmres without preconditioner');
% Nrestart_gmres = 10;
% Nit_gmres = 1000;
% [u_gmres, ~, relerr_gmres] = gmres(A, b, Nrestart_gmres, 1E-3, Nit_gmres);
% disp(['Relative residual gmres (source) ' num2str(relerr_gmres)]); %relerr_gmres = norm(A(u_gmres)-b) / norm(b);
% u_gmres = pagemtimes(sim.medium.Tr, reshape(u_gmres, size(u))); % remove the Tr^-1 scaling
% 
% %% BIGCSTAB
% % note: bicgstab does not work, operator A is called with
% % exponentially increasing inputs and the end result is 0!
% %
% disp('Running bicgstab without preconditioner');
% Nit_bicgstab = 100;
% [u_bicgstab, ~, relerr_bicgstab] = bicgstab(A, b, 1E-3, Nit_bicgstab);
% disp(['Relative residual bicgstab (source) ' num2str(relerr_bicgstab)]); %relerr_gmres = norm(A(u_gmres)-b) / norm(b);
% u_bicgstab = pagemtimes(sim.medium.Tr, reshape(u_bicgstab, size(u))); % remove the Tr^-1 scaling
% 
% %% Preconditioned
% Ap = @(x) sim.preconditioner(sim.operator(x));   % scaled operator L'+V'
% bp = sim.preconditioner(source.to_array());      % returns s' = Tl s
% 
% relerrs = norm(Ap(up(:))-bp) / norm(bp);
% disp(['Relative residual (source) ' num2str(relerrs)]);
% 
% %% Solve the problem using standard algorithms
% 
% %% GMRES
% disp('Running gmres with preconditioner');
% [u_gmresp, ~, relerr_gmresp] = gmres(Ap, bp, Nrestart_gmres, 1E-3, Nit_gmres);
% disp(['Relative residual gmres (source) ' num2str(relerr_gmresp)]); %relerr_gmres = norm(A(u_gmres)-b) / norm(b);
% u_gmresp = pagemtimes(sim.medium.Tr, reshape(u_gmresp, size(u))); % remove the Tr^-1 scaling
% 
% %% BIGCSTAB
% disp('Running bicgstab wit preconditioner');
% [u_bicgstabp, ~, relerr_bicgstabp] = bicgstab(Ap, bp, 1E-3, Nit_bicgstab);
% disp(['Relative residual bicgstab (source) ' num2str(relerr_bicgstabp)]); %relerr_gmres = norm(A(u_gmres)-b) / norm(b);
% u_bicgstabp = pagemtimes(sim.medium.Tr, reshape(u_bicgstabp, size(u))); % remove the Tr^-1 scaling
% 
% 
% %%
% figure;
% plot(z, u(4:4:end));
% hold on;
% plot(z, u_gmres(4:4:end));
% plot(z, u_gmresp(4:4:end));
% plot(z, u_bicgstab(4:4:end));
% plot(z, u_bicgstabp(4:4:end));
% 
% %plot(z, uprime(4:4:end));
% plot(z_inside + zstart, I_th, '--', 'LineWidth', 2.5);
% legend(sprintf('anysim %d iterations', state.iteration),...
%     sprintf('GMRES %d x %d iterations', Nrestart_gmres, Nit_gmres), ...
%     sprintf('GMRES %d x %d it + precond', Nrestart_gmres, Nit_gmres), ...
%     sprintf('BICGSTAB %d iterations', Nit_bicgstab), ...
%     sprintf('BICGSTAB %d it + precond', Nit_bicgstab), ...
%     'analytical');

