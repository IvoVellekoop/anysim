%   Simple test of the Pantograph toolbox
%   (c) 2022. Ivo Vellekoop
%   In this example, the Pantograph equation is not accretive and the method diverges.
%   However, by setting the option 'accretive=false', the Pantograph toolbox
%   solves the the anti-symmetrized system instead, which is accretive.
%
close all; clearvars;

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(10/opt.pixel_size);
opt.boundaries_width = 200;
opt.termination_condition = TerminationCondition(relative_limit= 1E-6);
opt.V_max = 0.5;
opt.callback = DisplayCallback();

%% Medium parameters
lambda = 0.9;
a = 0.1;%ones(opt.N,1);0.1;%(-0.1 + 2i) * ones(opt.N, 1); %(-0.1 + 2i) * ones(opt.N, 1);
b = -5;%-5*exp(0.1i * (1:opt.N));
%a(end-500:end) = 6; % 'manual' boundary conditions
%b(end-500:end) = 0;
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation. 
% First run the standard version, which will diverge because the
% equation is not accretive.
opt.accretive = true; % pretend that the equation is accretive, even though it is not
sim = PantographF(a, b, lambda, t0, opt);

% Define source
z = sim.grid.coordinates(1);
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z(1:t0));
source = sim.define_source(src);

%% Perform the different simulations and compare the results
[precond, ~] = compare_simulations(sim, source, []);%default_simulations);

%% Repeat with (anti-)symmetrized equation
s_opt = opt;
s_opt.accretive = false; % symmetrize system
s_sim = PantographF(a, b, lambda, t0, s_opt);
s_source = s_sim.define_source(src);
[s_precond, table] = compare_simulations(s_sim, s_source, []);%default_simulations);

%%
signal = [src; precond(end).value];
s_signal = [src; s_precond(end).value];
full_z = [z(1:t0); z+z(t0+1)];
%plot(full_z, real(signal));
%hold on;
plot(full_z, real(s_signal));
xlabel('t [s]');
ylabel('x(t)');
%legend('Anysim symmetrized', 'BiCGStab non-symmetrized');

%%
%inspect_sim(sim);

