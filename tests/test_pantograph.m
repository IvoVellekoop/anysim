%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(15/opt.pixel_size);
opt.boundaries_width = 0; % don't add boundaries
opt.callback = DisplayCallback();
opt = override_options(opt); % globally override options when calling this simulation from a script (see test_all)

%% Medium parameters
lambda = 1; %0.5; %1.5
a = 0.5;%+5i;
a = a * ones(1, opt.N);
a(100:200) = -0.3;
b = -0.1;%-4.9;
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
sim = PantographF(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
% Define source
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z(1:t0));
source = sim.define_source(src(:));

%% Perform the different simulations and compare the results
% note: the preconsitioned system has real eigenvalues, but it is
% non-normal and nonsymmetric.
simulations = default_simulations();
results = compare_simulations(sim, source, simulations);
%%
plot(z, results(end).value);
xlabel('t [s]');
ylabel('f(t)');

%%
%inspect_sim(sim);
