%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%
%   Simulates 1-D propagation of light through a slab of glass

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
opt.N = 256;
opt.boundaries_width = 64;
%opt.callback = DisplayCallback();
%opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
%opt.crop_to_roi = false;
%opt.legacy_mode = false; % makes both V and L accretive. No difference with new method if there is no absorption in the structure.


%% create an AnySim object for a homogeneous sample
n = ones(opt.N, 1); % refractive idex
n(100:130) = 1.5;
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(1); % intensity-only source (isotropic) at t=0

%% Compare to other methods and compute errors
simulations = default_simulations();

% without preconditioner, all methods diverge!
%bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations);




