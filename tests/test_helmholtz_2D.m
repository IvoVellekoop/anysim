%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Construct refractive index structure
% load the logo image, we will use the different color channels for
oversampling = 1;
im = imresize(single(imread("natlogo.png"))/255, oversampling, "bilinear");
n_iron = 2.8954 + 2.9179i;  % Fe at 532nm https://refractiveindex.info/?shelf=main&book=Fe&page=Johnson
n_gold = 0.54386 + 2.2309i;
n_material = n_gold;
%n = (im(:,:,3)-(im(:,:, 1)) > 0.25) * (n_gold-1) + 1;
n = (im(:,:,3)+(im(:,:, 1)) > 0.25) * (n_material-1) + 1;

%% Set up simulation options
opt = HelmholtzSimOptions();
tol = 1E-2;
opt.boundaries_width = 30; % 0=periodic
opt.callback = DisplayCallback();%plot_residual = true);
opt.termination_condition = TerminationCondition(relative_limit = tol);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*oversampling);

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(im(:,:,2));
[E, state] = sim.exec(source);

%%
I = abs(E).^2;
x0 = (1400-200) * sim.grid.pixel_size / 2;
close all; imagesc([-x0, x0], [-x0, x0], sqrt(I)); axis image; colorbar; colormap gray;
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
[precond, table] = compare_simulations(sim, source, simulations, tol = tol);

%% Compare with legacy method (wavesim)
l_opt = opt;
l_opt.legacy_mode = true;
l_sim = HelmholtzSim(n, l_opt);
l_source = l_sim.define_source(im(:,:,2));
[l_precond, l_table] = compare_simulations(l_sim, l_source, simulations, tol = tol);
table(end+1) = "legacy " + l_table(2);

