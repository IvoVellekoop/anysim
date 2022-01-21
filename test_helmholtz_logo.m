%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = struct(); % clear any previous options
im = imresize(single(imread("natlogo.png"))/255, 4);
n = 1.5 - 0.5 * im(:,:,1) + (0.5+0.4i) * im(:,:,3);
opt.N = [size(n,1), size(n,2), 1];   % Nx, Ny, Nz   (constant in z)
opt.boundaries.periodic = [false, false, true];
opt.boundaries.width = 64;
opt.callback.handle = @DisplayCallback;
opt.callback.show_boundaries = true;
opt.termination_condition.relative_limit = 1E-3;
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.pixel_size = 0.25/4;
opt.alpha = 0.75;

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(im(:,:,2));
[E, state] = sim.exec(source);

%%
close all; imagesc(log10(abs(E).^2), [-8 -2]); axis image; colorbar;
xlim([200 1400]);
ylim([200 1400]);
colormap hot;

%% Compare to other methods and compute errors
%% Perform the different simulations and compare the results
comp_opt.analytical_solution = 0;
comp_opt.tol = []; % []=stop when reached same residual as anysim
comp_opt.iter = 1000; %[]= never use more operator evaluations than anysim
simulations = default_simulations;

%comp_opt.preconditioned = false;
%bare = compare_simulations(sim, source, simulations, comp_opt);

%% Repeat with preconditioner
comp_opt.preconditioned = true;
%precond = compare_simulations(sim, source, simulations, comp_opt);

%[A, preA] = simulation_eigenvalues(sim);



