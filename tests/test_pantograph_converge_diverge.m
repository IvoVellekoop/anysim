%   Simple test of the Pantograph toolbox
%   (c) 2022. Ivo Vellekoop
%   In this example, the Pantograph equation is not accretive and the method diverges.
%   However, by setting the option 'accretive=false', the Pantograph toolbox
%   solves the the anti-symmetrized system instead, which is accretive.
%
close all; clear all;

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(10/opt.pixel_size);
opt.boundaries_width = 0; % don't add boundaries
opt.termination_condition = TerminationCondition(relative_limit= 1E-6);
opt.termination_condition_interval = 1;
opt.callback = DisplayCallback();
opt.V_max = 0.5;
%% Medium parameters
lambda = 0.5;
a = (0.1+5i) * ones(opt.N,1);
b = 5*cos(1:opt.N);
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
opt.accretive = false;true;
sim = PantographF(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
% Define source
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z(1:t0));
source = sim.define_source(src);

%% construct the analytical solution
analytical_solution = exp(-(a(1)+b(1)) .* z) * src(end);

%% Perform the different simulations and compare the results
%[comp, state] = sim.exec(source);
simulations = default_simulations();%has_adjoint=true);
[precond, table] = compare_simulations(sim, source, simulations, analytical_solution=analytical_solution);
comp = precond(end).value;%sim.exec(source);
%%
%plot(abs(comp));
plot(z, log(abs(comp)));
hold on;
%plot(abs(analytical_solution));
plot(z, log(abs(analytical_solution)));
hold off;
xlabel('t [s]');
ylabel('f(t)');
legend('computed', 'analytical');
%%
plot(comp(80:150));
hold on;
plot(analytical_solution(80:150));
%plot(abs(comp ./ analytical_solution));
%plot(log(abs(comp-analytical_solution)));


%%
%inspect_sim(sim);

