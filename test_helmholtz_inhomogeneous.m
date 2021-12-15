%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = struct(); % clear any previous options
opt.N = [256, 1, 1];   % Nx, Ny, Nz   (constant in z)
opt.boundaries.periodic = [false, true, true];
opt.boundaries.width = 64;
opt.callback.handle = @DisplayCallback;
opt.callback.show_boundaries = true;
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.pixel_size = 0.25;
opt.crop = false;

%% create an AnySim object for a homogeneous sample
n = ones(opt.N); % refractive idex
n(100:130) = 1.5;
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(ones(1,opt.N(2),1)); % intensity-only source (isotropic) at t=0
%[E, state] = sim.exec(source);


%% Compare to other methods and compute errors
%% Perform the different simulations and compare the results
comp_opt.analytical_solution = 0;
comp_opt.tol = []; % []=stop when reached same residual as anysim
comp_opt.iter = 1000; %[]= never use more operator evaluations than anysim
simulations = default_simulations;

comp_opt.preconditioned = false;
bare = compare_simulations(sim, source, simulations, comp_opt);

%% Repeat with preconditioner
comp_opt.preconditioned = true;
precond = compare_simulations(sim, source, simulations, comp_opt);

%[L, GL] = simulation_eigenvalues(sim);



