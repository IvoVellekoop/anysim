%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Construct refractive index structure
% load the logo image, we will use the different color channels for
oversampling = 1;
im = imresize(single(imread("natlogo.png"))/255, 0.25, "bilinear");
n_iron = 2.8954 + 2.9179i;  % Fe at 532nm https://refractiveindex.info/?shelf=main&book=Fe&page=Johnson
n_gold = 0.54386 + 2.2309i;
n_material = n_iron;

%n = (im(:,:,3)+(im(:,:, 1)) > 0.25) * (n_material-1) + 1;

%% Set up simulation options
opt = HelmholtzSimOptions();
opt.boundaries_width = 20; % 0=periodic
%opt.callback = DisplayCallback();%plot_residual = true);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.termination_condition = TerminationCondition(relative_limit = 1E-5);
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*oversampling);
opt.N = [100, 100];

n = zeros(opt.N);
grid_center = opt.N/2;
positions = [0 0; 1 1; -1 -1] * 20 + grid_center;
radius = 10;
x = 1:opt.N(1);
y = (1:opt.N(2)).';
for pos_idx = 1:size(positions, 1)
    n = n + ((x - positions(pos_idx, 1)).^2 + (y - positions(pos_idx, 2)).^2 < radius^2);
end
n = 1 + (n_iron - 1) * n;
n_material = n_gold;
n = (im(:,:,3)-(im(:,:, 1)) > 0.25) * (n_material-1) + 1;

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
src = exp(-((-3:3).^2 + (-3:3)'.^2)); 
source = sim.define_source(src, grid_center-3);
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
simulations = default_simulations("nonsymmetric");

% without preconditioner, all methods diverge!
%bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations, tol = tol);

%% Compare with legacy method (wavesim)
l_opt = opt;
l_opt.legacy_mode = true;
l_sim = HelmholtzSim(n, l_opt);
l_source = l_sim.define_source(src, grid_center-3);
[l_precond, l_table] = compare_simulations(l_sim, l_source, simulations, tol = tol);
table(end+1) = "legacy " + l_table(2);

%%
for i=1:length(precond)
    fprintf("%f %f\n", norm(precond(i).value(:) - E(:)), norm(l_precond(i).value(:) - E(:)));
end
