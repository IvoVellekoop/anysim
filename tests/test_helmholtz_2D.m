%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Construct refractive index structure
% load the logo image, we will use the different color channels for
oversampling = 1;
im = imresize(single(imread("natlogo.png"))/255, oversampling, "bilinear");
n_iron = 2.8954 + 2.9179i;  % Fe at 532nm https://refractiveindex.info/?shelf=main&book=Fe&page=Johnson
n = (mean(im, 3) > 0.25) * n_iron;

%% Set up simulation options
opt = HelmholtzSimOptions();
opt.boundaries_width = 0; %periodic
opt.callback = DisplayCallback();
%opt.termination_condition = TerminationCondition(relative_limit = 1E-3);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*oversampling);
opt.alpha = 0.75;

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
% source = sim.define_source(im(:,:,2));
source = im(:,:,1) + 1i * im(:,:,2);
[E, state] = sim.exec(source);

%%
I = abs(E).^2;
x0 = (1400-200) * sim.grid.pixel_size / 2;
close all; imagesc([-x0, x0], [-x0, x0], sqrt(I), [0 0.4E-0]); axis image; colorbar;
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

%% Compare with legacy method (wavesim)
l_opt = opt;
l_opt.legacy_mode = true;
l_sim = HelmholtzSim(n, l_opt);
l_source = l_sim.define_source(im(:,:,2));
[l_precond, l_table] = compare_simulations(l_sim, l_source, simulations);
table(end+1) = "legacy " + l_table(2);

