%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%   Test case for lambda = 1, where we have an analytical solution

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(10/opt.pixel_size);
opt.boundaries_width = 0; % don't add boundaries
opt.termination_condition = TerminationCondition(relative_limit= 1E-6);
opt.V_max = 0.5;
%% Medium parameters
lambda = 1;
a = 3.5 + 1i;
b = 1 + 1i; % note, factor sqrt(lambda) included in beta to make Λ unitary
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
sim = Pantograph(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
% Define source
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z(1:t0-1));
source = sim.define_source(src(:));

%% Perform the different simulations and compare the results
zdil = z(t0:end);
zdil = zdil - zdil(1) + opt.pixel_size * 0.5;
analytical_solution = zeros([opt.N, 1]);
analytical_solution(1:t0-1) = src;
analytical_solution(t0:end) = exp(-(a+b) * zdil) * src(end);

simulations = default_simulations("nonsymmetric");
[precond, table] = compare_simulations(sim, source, simulations, analytical_solution=analytical_solution);
comp = precond(end).value / sim.Tr;%sim.exec(source);
%%
%plot(abs(comp));
plot(log(abs(comp)));
hold on;
%plot(abs(analytical_solution));
plot(log(abs(analytical_solution)));
hold off;
xlabel('t [s]');
ylabel('f(t)');
legend('computed', 'analytical');
%%
plot(abs(comp ./ analytical_solution));
%plot(log(abs(comp-analytical_solution)));

%%
%[~, GL] = simulation_eigenvalues(sim);

