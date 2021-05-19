%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%

% Note: in this simulation, all axes (including the time axis)
% are in micrometers. The absorbtion coefficient and the diffusion
% coeffecient both given in unit [um] and [1/um], respectively.
% To scale to time in seconds and diffusion coefficient in um^2/s just
% divide or multiply by c (in um/s)
%

% Set up size of simulation domain and number of grid points in x,y,z,t 
% dimensions.
opt = struct(); % clear any previous options
opt.pixel_size = {1 'um'};
opt.N = [256*4, 1, 1, 1]; %Nx, Ny, Nz, t   (constant in z and t)
opt.boundaries.periodic = true; %we manually define the boundaries inside the simulation domain
opt.callback.handle = @DisplayCallback;
opt.callback.cross_section = @(u) u(4,:,:);
opt.termination_condition.relative_limit = 1E-9;

%% Construct medium 
N_boundary = 200;
Dslab = 2;
a = zeros(opt.N(1), 1);                 % absorption coefficient (1/m)
D = ones(opt.N(1), 1, 1, 1) * Dslab;    % diffusion length (m)


%%
sz = size(D);
R = 0.5;
ze = 2*(1+R)/(1-R)*Dslab;
a(1:N_boundary) = Dslab/ze^2;  % boundary conditions
a((end-N_boundary+1):end) = Dslab/ze^2;
%D(1:zl-1) = Inf;
%D(zr+1:end) = Inf;
%a(1:zl-1) = 0;
%a(zr+1:end) = 0;

%%
sim = DiffuseSim(D, a, opt);
z = sim.grid.crop(sim.grid.coordinates(1));

%% Define source and run the simulation
ell = Dslab*3;
z_indices = N_boundary+1:opt.N(1)-N_boundary;
z_inside = z(z_indices);
zl = z_inside(1) - 0.5 * opt.pixel_size{1};
zr = z_inside(end) + 0.5 * opt.pixel_size{1};
zstart = z_inside(1);
z_inside = z_inside - zstart;
Isource = exp(-(z_inside.'+ 0.5 * opt.pixel_size{1})/ell)/ell;
source = sim.define_source(Isource, [4, N_boundary+1,ceil(opt.N(2)/2),1,1]); % intensity-only source (isotropic) at t=0

u = sim.exec(source);

%%
Iz = squeeze(u(4,:,ceil(end/2)));
Fz = squeeze(u(1,:,ceil(end/2)));
figure;
plot(z, Iz, 'LineWidth', 1.5);
hold on;

% theoretical intensity distribution (normalized to same height as simulation):
L = zr-zl;
I0 = 1/Dslab/(L+2*ze);
e1 = exp(-z_inside/ell);
e2 = exp(-L/ell);
I_th = I0 * (-e1 * (L + 2*ze) * ell + e2 * (z_inside + ze) * (ell - ze) + (L - z_inside + ze) * (ze + ell));

ylimits = [0, 10];
plot(z_inside + zstart, I_th, '--', 'LineWidth', 2.5);
plot([zl, zl], ylimits, 'k:', 'LineWidth', 1); 
plot([zr, zr], ylimits, 'k:', 'LineWidth', 1); 
ylim(ylimits);
xlabel('z [\mum]');
ylabel('u [J/m^3]');
set(gca, 'FontSize', 14);
legend('simulated', 'theoretical', 'sample boundaries', 'FontSize', 12);
title('Intensity profile inside diffusing slab')
%
relerr = norm(Iz(z_indices).'-I_th(:))/norm(I_th(:));
disp(['Relative error ' num2str(relerr)]);
