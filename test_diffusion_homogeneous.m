%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%
% This test simulates steady state diffusion in a 1-D homogeneous medium
% with absorption, and compares the result to the analytical solution
%

% Set up size of simulation domain and number of grid points in x,y,z,t 
% dimensions.
N = 256;
opt.pixel_size = {0.5 'um'};
opt.N = [N, 1, 1, 1];
opt.boundaries.periodic = true; % all boundaries periodic
opt.gpu_enabled = false;
opt.termination_condition.handle = @tc_fixed_iteration_count;
opt.termination_condition.iteration_count = 100;
opt.callback.interval = 10;
opt.callback.handle = @DisplayCallback;
opt.callback.cross_section = @(u) u(4,:);
opt.V_max = 0.95;

%% Construct medium 
a = 0.1; % absorption coefficient (1/um)
D = 2; % diffusion coefficient

sim = DiffuseSim(D, a, opt);

%% Place a source in center and run simulation
source = sim.define_source(ones(1,1,opt.N(2)), [4,ceil(N/2),1,1,1]);
u = sim.exec(source);
I = shiftdim(u(4,:), 1);

%% Compute theoretical solution
x = sim.coordinates(1);
x = x - (ceil(N/2)-1) * (x(2)-x(1)); % shift coordinates relative to source
semilogy(x, I)
hold on;

mueff = sqrt(a/D);
h = opt.pixel_size{1,1};
hs = 1.0i * pi/h;

% Analytical solution with sinc source.
xa = abs(x);
xa(xa<1E-6)=1E-6;
I_th = h*mueff / (4*pi*a) * (...
    2*pi*exp(-mueff*xa)...
    + 1.0i*exp(-mueff*xa) .* (expint(xa*(-hs-mueff)) - expint(xa*(hs-mueff)))...% correction for finite pixel size
    + 1.0i*exp(mueff*xa) .* (-expint(xa*(-hs+mueff)) + expint(xa*(hs+mueff))));
% note, the expression above differs from the one in mathematica because 
% of the different definition of expint. The sign of x is reversed, as is
% the sign of the expint itself. Moreover, i pi is added or subtracted 
semilogy(x, I_th); 

%% Determine accuracy
dI = I - I_th;
relative_error = norm(dI(:))^2/norm(I(:))^2;
disp(relative_error)

