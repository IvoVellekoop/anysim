%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
opt.N = 256;
opt.boundaries_width = 64;
opt.callback = DisplayCallback();
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.crop_to_roi = false;
opt.alpha = 0.75;

%% create an AnySim object for a homogeneous sample
n = ones(opt.N, 1); % refractive idex
n(100:130) = 1.5;
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(1); % intensity-only source (isotropic) at t=0
%[E, state] = sim.exec(source);


%% Compare to other methods and compute errors
%% Perform the different simulations and compare the results
simulations = default_simulations;
bare = compare_simulations(sim, source, simulations, preconditioned=false);

%% Repeat with preconditioner
precond = compare_simulations(sim, source, simulations);

%[A, preA] = simulation_eigenvalues(sim);



