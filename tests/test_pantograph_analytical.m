%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%   Test case for lambda = 1, where we have an analytical solution
close all; clear all;

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(10/opt.pixel_size);
opt.gpu_enabled = false;
opt.boundaries_width = 100; % absorbing boundaries
opt.termination_condition = TerminationCondition(relative_limit= 1E-6);
opt.termination_condition_interval = 1;
opt.callback = DisplayCallback();
opt.V_max = 0.5;
%% Medium parameters
lambda = 1;
a = 2 * ones(opt.N, 1);
b = -1i * ones(opt.N, 1);
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
opt.accretive = true;
sim = PantographF(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
% Define source
src = zeros(t0, 1);
src(end) = 1;
source = sim.define_source(src);

%% construct the analytical solution
analytical_solution = exp(-(a(1)+b(1)) .* z) * src(end);

%% Perform the different simulations and compare the results
%[comp, state] = sim.exec(source);
simulations = default_simulations();
[precond, table] = compare_simulations(sim, source, simulations, analytical_solution=analytical_solution);
comp = precond(end).value;%sim.exec(source);
%%
plot(z, log(abs(comp)));
hold on;
%plot(abs(analytical_solution));
plot(z, log(abs(analytical_solution)));
hold off;
xlabel('t [s]');
ylabel('f(t)');
legend('computed', 'analytical');

%%
%inspect_sim(sim);

