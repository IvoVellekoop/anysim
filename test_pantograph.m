%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Simulation parameters
opt = struct(); % clear any previous options
opt.pixel_size = 0.01;
opt.pixel_unit = 's';
opt.N = [256*4, 1, 1, 1]; %Nx, Ny, Nz, t   (constant in z and t)
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.dilation_start = 100; % first values are starting condition
opt.callback.handle = @DisplayCallback;
opt.boundaries.periodic = true; % don't add boundaries
opt.gpu_enabled = false;

%% Medium parameters
s = 0.5;
a = -10;
%b = zeros(opt.N);
b = -5; %0

%% Set up AnySim simulation
sim = Pantograph(a, b, s, opt);
z = sim.coordinates(1);
zdil = z(opt.dilation_start:end);
zdil = zdil - zdil(1) + opt.pixel_size;
% Define source
src = 0.1+exp(-0.5*((z-0.85)/.02).^2)-0.5*exp(-0.5*((z-0.80)/.05).^2);
src = src(1:opt.dilation_start-1);
%src = 1:opt.dilation_start-1;
source = sim.define_source(src(:));

%% Perform the different simulations and compare the results
comp_opt.analytical_solution = zeros([opt.N, 1]);
comp_opt.analytical_solution(1:opt.dilation_start-1) = src;
comp_opt.analytical_solution(opt.dilation_start:end) = exp((a+b) * zdil) * src(end);
comp_opt.tol = [];

simulations = default_simulations;
comp_opt.preconditioned = true;
precond = compare_simulations(sim, source, simulations, comp_opt);

%[L, GL] = simulation_eigenvalues(sim);
