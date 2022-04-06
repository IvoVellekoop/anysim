%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
oversampling = 1;
n_iron = 2.8954 + 2.9179i;  % Fe at 532nm https://refractiveindex.info/?shelf=main&book=Fe&page=Johnson
%im = imresize(single(imread("natlogo.png"))/255, oversampling, "bilinear");
%n = 1 + n_iron * im(:,:,3) * 0.2; % iron
y = sim.grid.coordinates(1);
x = sim.grid.coordinates(2);
grid_extent = sim.grid.pixel_size .* sim.grid.N;
grid_center = [x(1+floor(end/2)), y(1+floor(end/2))];
im = exp(-0.5 * ((((x - grid_center(1)) - -0.45 * grid_extent(1)) ./ 1.0).^2 + ((y - grid_center(2)) ./ 10).^2));  % The source position
im = im .* exp(2i * pi * x);
im = repmat(im, [1, 1, 3]);
im(:,:,1) = real(im(:,:,1));
im(:,:,2) = imag(im(:,:,2));
im(:,:,3) = 0;  % The object
positions = [0 0; 1 1; -1 -1] * 10 + grid_center;
radius = 4;
for pos_idx = 1:size(positions, 1)
    im(:,:,3) = im(:,:,3) + ((x - positions(pos_idx, 1)).^2 + (y - positions(pos_idx, 2)).^2 < radius^2);
end
n = 1 + (n_iron - 1) * im(:,:,3) * 0.2;  % A fraction of that of iron
boundary_width = 0.1;
n = n + 0.2i * max(0, max(abs(x - grid_center(1)) / grid_extent(1), abs(y - grid_center(2)) / grid_extent(2)) - (0.5 - boundary_width)) ./ boundary_width;  % Add absorbing boundary
opt.N = [size(n,1), size(n,2)];   % Nx, Ny, Nz   (constant in z)
opt.boundaries_width = 0; %periodic (Setting this to > 0 requires nuttallwin function from Signal Processing Toolbox)
% opt.callback = DisplayCallback();
opt.termination_condition = TerminationCondition(relative_limit = 1E-3);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.pixel_size = 0.25/oversampling;
opt.alpha = 0.75;
opt.crop_to_roi = false; % so that we can compare with the forward operator

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

