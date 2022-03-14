%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(5/opt.pixel_size);
%opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.callback = DisplayCallback();
opt.boundaries_width = 0; % don't add boundaries

%% Medium parameters
lambda = 0.5; %1.5
a = 10;
b = 5 / sqrt(lambda); % note, factor sqrt(lambda) included in beta to make Λ unitary
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
sim = Pantograph(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
zdil = z(t0:end);
zdil = zdil - zdil(1) + opt.pixel_size;
% Define source
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z);
src = src(1:t0-1);
%src = 1:opt.dilation_start-1;
source = sim.define_source(src(:));

%% Perform the different simulations and compare the results
comp_opt.analytical_solution = zeros([opt.N, 1]);
comp_opt.analytical_solution(1:t0-1) = src;
comp_opt.analytical_solution(t0:end) = exp(-a * zdil) * src(end);
comp_opt.tol = [];

% note: the preconsitioned system has real eigenvalues, but it is
% non-normal and nonsymmetric.
simulations = default_simulations("nonsymmetric");
%no_precond = compare_simulations(sim, source, simulations, preconditioned = false);
precond = compare_simulations(sim, source, simulations);

%%
[L, GL] = simulation_eigenvalues(sim);
