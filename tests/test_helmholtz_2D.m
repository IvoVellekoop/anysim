%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%
%   Simulates propagation of light in a 2-D structure made of iron (uses the scalar wave equation)

%% Construct refractive index structure
% load the logo image
oversampling = 0.25;
im = single(imread("logo_structure_vector.png"))/255;
n_gold = 0.54386 + 2.2309i; % Au at 532nm https://refractiveindex.info/?shelf=main&book=Au&page=Johnson
n_silver = 0.054007 + 3.4290i;
n_copper = 1.1159 + 2.5956i;
n_iron = 2.8954 + 2.9179i;
n_contrast = n_iron - 1;

% this version has a slight dominance for the new method, at Rich60.
%n = (im(:,:,3)-(im(:,:, 1)) > 0.25) * n_contrast + 1;

% this version has a very clear dominance for the new method, at Rich80.
% difference may be caused by treatment of boundaries?
n = imresize((im(:,:,3) > 0.25) * n_contrast + 1, oversampling, "bilinear");
src = imresize(im(:,:,2), oversampling, "bilinear");

%% Set up simulation options
opt = HelmholtzSimOptions();
%tol = 0.005;
opt.boundaries_width = 30; % 0=periodic
%opt.callback = DisplayCallback();%plot_residual = true);
opt.termination_condition = TerminationCondition(relative_limit = 1E-4, iteration_count = 1E6);
opt.wavelength = 0.532;
opt.pixel_size = opt.wavelength/(3*max(abs(n_contrast+1)));
opt.gpu_enabled = true;

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
%h = gcf();
%h.PaperUnits = 'centimeters';
%h.PaperSize = [1, 1];
%h.PaperPositionMode = 'manual';
xticks([]);
yticks([]);
axis on;
print('helmholtz_logo.pdf', '-dpdf', '-fillpage'); %cannot use eps because of the alpha channel

%% Compare to other methods and compute errors
simulations = default_simulations();

% without preconditioner, all methods diverge!
%bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations, iter=1E5, measure_time=false);

%% Compare with legacy method (wavesim)
l_opt = opt;
l_opt.legacy_mode = true;
l_sim = HelmholtzSim(n, l_opt);
l_source = l_sim.define_source(src);
[l_precond, l_table] = compare_simulations(l_sim, l_source, simulations, iter=1E5);
table(end+1) = "legacy " + l_table(2);

%% 
fprintf("Relative error compared to accurate simulation\n");
nE = norm(E(:));
for s = 1:length(precond)-1
    fprintf("%s %f %f\n", precond(s).name, norm(precond(s).value(:)-E(:))/nE, norm(l_precond(s).value(:)-E(:))/nE);
end
