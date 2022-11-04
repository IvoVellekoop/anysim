%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%
%   Simulates propagation of light in a 2-D structure made of iron (uses the scalar wave equation)

%% Construct refractive index structure
% load the logo image
oversampling = 0.25;
im = single(imread("logo_structure_vector.png"))/255;
n_water = 1.33;
n_fat = 1.46;

n = imresize((im(:,:,3) > 0.25) * (n_fat-n_water) + n_water, oversampling, "bilinear");
src = imresize(im(:,:,2), oversampling, "bilinear");

%% Set up simulation options
opt = HelmholtzSimOptions();
opt.boundaries_width = 75;
opt.termination_condition = TerminationCondition(tolerance = 1E-3, iteration_count = 1E4);
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*abs(n_fat));
opt.gpu_enabled = true;
opt.callback = DisplayCallback();
opt = override_options(opt); % globally override options when calling this simulation from a script (see test_all)

%% create the AnySim object
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(src);
[E, state] = sim.exec(source);

%%
I = abs(E);
x = sim.grid.coordinates(1); x = x - x(end)/2;
close all; 
figure;
imshow(I/max(I(:)), [], Colormap=parula(256), XData = x, YData = x); axis image;
cb = colorbar; cb.FontSize = 18;
colormap(parula);
hold on;
mask = imgaussfilt(single(edge(im(:,:,3), 'canny')), 2); mask = sqrt(mask) / max(mask(:));
mask(1:10,:) = 0; mask(:,1:10) = 0;
mask(end-10:end,:) = 0; mask(:, end-10:end) = 0;
overlay = imshow(ones([size(mask), 3]), [], XData = x, YData = x);
overlay.AlphaData = mask;
hold off;
scalebar_width = 2.5; % in micrometer
rectangle("Position", [x(round(length(x) * 0.85)), x(round(length(x) * 0.95)), scalebar_width, scalebar_width/5], FaceColor='w', LineWidth=0.5);
xticks([]);
yticks([]);
axis on;

%% Compare to other methods and compute errors
simulations = default_simulations();
results = compare_simulations(sim, source, simulations);


