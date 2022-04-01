%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
oversampling = 2;
im = imresize(single(imread("natlogo.png"))/255, oversampling, "bilinear");
%n = 1.5 - 0.5 * im(:,:,1);% + (0.5+0.4i) * im(:,:,3);
n = 1 + (1.8 + 2.9i) * im(:,:,3) * 0.2; % iron
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
I = abs(E).^2;
x0 = (1400-200) * sim.grid.pixel_size / 2;
close all; imagesc([-x0, x0], [-x0, x0], sqrt(I), [0 0.4E-1]); axis image; colorbar;
xlabel('x [\mum]');
ylabel('y [\mum]');
print(gcf, '-depsc', 'helmholtz_logo.eps');
%xlim([-25 25]);
%ylim([-25 25]);
%colormap hot;

%% Compare to other methods and compute errors
simulations = default_simulations("nonsymmetric", has_adjoint = true);

% without preconditioner, all methods diverge!
%bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations);

