%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%

% Note: in this simulation, all axes (including the time axis)
% are in micrometers. The absorbtion coefficient and the diffusion
% coeffecient both given in unit [um] and [1/um], respectively.
% To scale to time in seconds and diffusion coefficient in um^2/s just
% divide or multiply by c (in um/s)
%

%% Simulation parameters
opt = DiffuseSimOptions();
opt.N = 1024; %1-D simulation
opt.pixel_size = 1;
opt.pixel_unit = 'um';
opt.boundaries_width = 0; %we manually define the boundaries inside the simulation domain
opt.callback = DisplayCallback(component = 4);
opt.termination_condition = TerminationCondition(relative_limit = 1E-5);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms


%% Medium parameters
N_boundary = 200;   % boundary width in pixels
Dslab = 2;          % diffusion length (m)
aslab = 0;          % absorption coefficient (1/m)
R = 0.5;            % angle-averaged reflection coefficient at interface.

%% Construct medium 
D = ones(opt.N, 1, 1, 1) * Dslab;
a = zeros(opt.N, 1) * aslab;

% construct boundary outside medium, with absorption coefficient
% tuned to match mixed boundary conditions with extrapolation length
%
ze = 2*(1+R)/(1-R)*Dslab;       % extrapolation length
a(1:N_boundary) = Dslab/ze^2;   % fill boundaries with matching absorption
a((end-N_boundary+1):end) = Dslab/ze^2;

%% Set up AnySim simulation with an exponentially decaying source
sim = DiffuseSim(D, a, opt);

% Define source
z = sim.grid.coordinates(1);

ell = Dslab*3;
z_indices = N_boundary+1:opt.N-N_boundary;   % indices of slab
z_inside = z(z_indices);                        % corresponding z-coordinates
zl = z_inside(1) - 0.5 * opt.pixel_size;        % z-coordinate of left boundary
zr = z_inside(end) + 0.5 * opt.pixel_size;      % z-coordinate of right boundary
zz = z_inside-zl;                               % position with respect to sample start
Isource = exp(-zz'/ell)/ell;                    % exponentially decaying source
source = sim.define_source(Isource, [4, N_boundary+1]); % intensity-only source (isotropic) at t=0

%% theoretical intensity distribution (normalized to same height as simulation):
L = zr-zl;
I0 = 1/Dslab/(L+2*ze);
e1 = exp(-zz/ell);
e2 = exp(-L/ell);
I_th = I0 * (-e1 * (L + 2*ze) * ell + e2 * (zz + ze) * (ell - ze) + (L - zz + ze) * (ze + ell));


%% Perform the different simulations and compare the results
analytical_solution = nan([4, opt.N]);
analytical_solution(4, z_indices) = I_th(:);

%% Compare to other methods and compute errors
simulations = default_simulations(has_adjoint = true);

% without preconditioner, all methods diverge!
bare = compare_simulations(sim, source, simulations, preconditioned = false, analytical_solution=analytical_solution);

%%
[precond, table] = compare_simulations(sim, source, simulations, analytical_solution=analytical_solution);
%%
[L, GL] = simulation_eigenvalues(sim);

%%
semilogy(z, precond(end).value(4,:));
hold on;
semilogy(z, analytical_solution);
hold off;

