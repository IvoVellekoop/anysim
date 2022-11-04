%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%
%   Simulates 1-D propagation of light through a slab of glass

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
opt.N = 256;
opt.boundaries_width = 64;
opt.forward_operator = true;
opt.preconditioner = "moborn";%"shift";
opt.preconditioner_shift = 1;
%opt.callback = DisplayCallback();
%opt.callback_interval = 100;
opt = override_options(opt); % globally override options when calling this simulation from a script (see test_all)


%% create an AnySim object for a homogeneous sample
n = ones(opt.N, 1); % refractive idex
n(100:130) = 1.5;
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(1); % intensity-only source (isotropic) at t=0

%% Compare to other methods and compute errors
simulations = default_simulations();
results = compare_simulations(sim, source, simulations);




