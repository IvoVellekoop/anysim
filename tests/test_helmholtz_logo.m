%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
oversampling = 4;
im = imresize(single(imread("natlogo.png"))/255, oversampling, "bilinear");
%n = 1.5 - 0.5 * im(:,:,1);% + (0.5+0.4i) * im(:,:,3);
n = 1 + (1.8 + 2.9i) * im(:,:,1) * 0.2; % iron
opt.N = [size(n,1), size(n,2)];   % Nx, Ny, Nz   (constant in z)
opt.boundaries_width = 0; %periodic
opt.callback = DisplayCallback();
opt.termination_condition = TerminationCondition(relative_limit = 1E-3);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.pixel_size = 0.25/oversampling;
opt.alpha = 0.75;
opt.crop_to_roi = false; % so that we can compare with the forward operator

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(im(:,:,2));
[E, state] = sim.exec(source);

%%
I = abs(E(264:1464, 264:1464)).^2;
x0 = (1400-200) * sim.grid.pixel_size / 2;
close all; imagesc([-x0, x0], [-x0, x0], log10(I), [-8 -2]); axis image; colorbar;
xlabel('x [\mu m]')
ylabel('y [\mu m]')
colormap hot;

%% Compare to other methods and compute errors
simulations = default_simulations("nonsymmetric");

% without preconditioner, all methods diverge!
%bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations);

