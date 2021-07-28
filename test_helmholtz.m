%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = struct(); % clear any previous options
opt.N = [256, 1, 1];   % Nx, Ny, Nz   (constant in z)
opt.boundaries.periodic = [false, true, true];
opt.boundaries.width = 2*512;
opt.callback.handle = @DisplayCallback;
opt.callback.show_boundaries = true;
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.pixel_size = 0.25;
opt.crop = false;
opt.V_max = 0.618034; % optimal in limit infitely large simulation: norm(E)->inf
%0.562 (optimal for current situation: norm(E)/norm(source)=0.76;
%opt.V_max = 0.95;

%% create an AnySim object for a homogeneous sample
n = 0.9; % refractive idex
sim = HelmholtzSim(n, opt);

%% Define source and run the simulation
source = sim.define_source(ones(1,opt.N(2),1)); % intensity-only source (isotropic) at t=0
[E, state] = sim.exec(source);

%% calculate exact solution analytically
k = n*2*pi/sim.opt.wavelength;
x = abs(sim.coordinates(1));
h = sim.grid.pixel_size(1);
% To determine peak value at x=0, realize that Ei(x) ~ ln(x) for x close to
% 0 and find E(0) = -h/(4*pi*k)*2*log((pi/h+k)/(pi/h-k))  + i*h/(2*k)
% which can be rewritten as E(0) = -h/(2*pi*k)*log((1+h*k/pi)/(1-h*k/pi)) + i*h/(2*k)
% E(0) = -h/(pi*k) * atanh(h*k/pi) + i*h/(2*k)
% E(0) = h/(pi*k) * (i*pi/2 -atanh(h*k/pi))
% E(0) = i*h/(2*k) * (1 + 2i*atanh(h*k/pi)/pi)
% note: h*k/pi = 2*h/lambda  ==> Sampling rate with respect to Nyquist
% note: the Ei part (sol-plane wave component) is real and rapidly
% oscillating
% (note: bug in wavesim takes the other solution (where absorption is in
% the positive imaginary part)
    x = abs(x);
    phi = k * x;
    E_theory = 1.0i*h/(2*k)*exp(1.0i * phi)... %<--propagating plane wave.
        -h/(4*pi*k) * (...
        exp(1.0i * phi) .* (  expint(1.0i * (k-pi/h) * x) - expint(1.0i * (k+pi/h) * x)) -...
        exp(-1.0i * phi) .* ( -expint(-1.0i * (k-pi/h)* x)  + expint(-1.0i* (k+pi/h) * x)));

    % special case for values close to 0
    small = abs(k*x)<1E-10;
    E_theory(small) = 1.0i * h/(2*k) * (1+2i*atanh(h*k/pi)/pi); %exact value at 0.

%% Compare to other methods and compute errors
%% Perform the different simulations and compare the results
comp_opt.analytical_solution = sim.grid.pad(E_theory, 0, nan);
comp_opt.tol = []; % []=stop when reached same residual as anysim
comp_opt.iter = 1000; %[]= never use more operator evaluations than anysim
simulations = default_simulations;

comp_opt.preconditioned = false;
bare = compare_simulations(sim, source, simulations, comp_opt);

%% Repeat with preconditioner
comp_opt.preconditioned = true;
precond = compare_simulations(sim, source, simulations, comp_opt);





