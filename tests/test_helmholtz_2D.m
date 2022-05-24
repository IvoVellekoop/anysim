%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Construct refractive index structure
% load the logo image, we will use the different color channels for
oversampling = 1;
im = imresize(single(imread("natlogo.png"))/255, 0.75*oversampling, "bilinear");
n_gold = 0.54386 + 2.2309i; % Au at 532nm https://refractiveindex.info/?shelf=main&book=Au&page=Johnson
n = (im(:,:,3)-(im(:,:, 1)) > 0.25) * (n_gold-1) + 1;
%n = (im(:,:,3) > 0.25) * (n_gold-1) + 1;
%n(im(:,:,1) > 0.2) = 1 + 0.2i;


%% Set up simulation options
opt = HelmholtzSimOptions();
tol = 0.005;
opt.boundaries_width = 30; % 0=periodic
%opt.callback = DisplayCallback();%plot_residual = true);
opt.termination_condition = TerminationCondition(relative_limit = tol / 10, iteration_count = 1E6);
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*oversampling);

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(im(:,:,2));
[E, state] = sim.exec(source);

%%
I = abs(E).^2;
x = sim.grid.coordinates(1);
close all; imagesc(x-x(end)/2, x-x(end)/2, I.^0.25); axis image; colorbar; colormap gray;
xlabel('x [\mum]');
ylabel('y [\mum]');
print(gcf, '-depsc', 'helmholtz_logo.eps');
%xlim([-25 25]);
%ylim([-25 25]);
colormap(viridis());

%% Compare to other methods and compute errors
simulations = default_simulations("nonsymmetric");

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

%% 
fprintf("Relative error compared to accurate simulation\n");
nE = norm(E(:));
for s = 1:length(precond)-1
    fprintf("%s %f %f\n", precond(s).name, norm(precond(s).value(:)-E(:))/nE, norm(l_precond(s).value(:)-E(:))/nE);
end
