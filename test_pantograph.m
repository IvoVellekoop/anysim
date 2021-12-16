%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Simulation parameters
opt = struct(); % clear any previous options
opt.pixel_size = 0.001;
opt.pixel_unit = 's';
opt.N = [256*4, 1, 1, 1]; %Nx, Ny, Nz, t   (constant in z and t)
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.dilation_start = 17; % first 16 values are starting condition
opt.callback.handle = @DisplayCallback;
opt.boundaries.periodic = true; % don't add boundaries
opt.gpu_enabled = false;

%% Medium parameters
s = 0.5;
a = -10;
b = 0;%5;
src = 1:opt.dilation_start-1;

%% Set up AnySim simulation with an exponentially decaying source
sim = Pantograph(a, b, s, opt);
z = sim.coordinates(1);
zdil = z(opt.dilation_start:end);
zdil = zdil - zdil(1) + opt.pixel_size;
% Define source
source = sim.define_source(src);

%% Perform the different simulations and compare the results
comp_opt.analytical_solution = zeros([1, opt.N]);
comp_opt.analytical_solution(1, 1:opt.dilation_start-1) = src;
comp_opt.analytical_solution(opt.dilation_start:end) = exp(a * zdil) * src(end);
comp_opt.tol = [];

simulations = default_simulations;
comp_opt.preconditioned = true;
precond = compare_simulations(sim, source, simulations, comp_opt);

%[L, GL] = simulation_eigenvalues(sim);
