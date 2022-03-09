%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%
% This test simulates steady state diffusion in a 1-D homogeneous medium
% with absorption, and compares the result to the analytical solution
%

%% Set up simulation options
opt = DiffuseSimOptions();
opt.N = 256;
opt.grid.boundaries_width = 0; % all boundaries periodic
opt.grid.pixel_size = 0.5;
opt.grid.pixel_unit = 'um';
opt.callback = DisplayCallback(cross_section = {4});
opt.termination_condition_interval = 1;
opt.forward_operator = true;

%% Construct medium 
a = 0.1;    % absorption coefficient [um^-1]
D = 2;      % diffusion coefficient [um]

sim = DiffuseSim(D, a, opt);

%% Place a source in center and run simulation
source = sim.define_source(ones(1,1,1), [4,ceil(opt.N(1)/2),1,1,1]);
[u, state] = sim.exec(source);
I = shiftdim(u(4,:), 1);

%% Compute theoretical solution
x = sim.grid.coordinates(1);
x = x - (ceil(opt.N(1)/2)-1) * (x(2)-x(1)); % shift coordinates relative to source
semilogy(x, I)
hold on;

mueff = sqrt(a/D);
h = opt.grid.pixel_size;
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
u_th = zeros(4, length(dI)) + nan;
u_th(4, :) = I_th;
simulations = default_simulations;
bare = compare_simulations(sim, source, simulations, preconditioned = false, analytical_solution = u_th);

%% Repeat with preconditioner
precond = compare_simulations(sim, source, simulations, analytical_solution = u_th);



