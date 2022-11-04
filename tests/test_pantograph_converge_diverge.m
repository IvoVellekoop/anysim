%   Simple test of the Pantograph toolbox
%   (c) 2022. Ivo Vellekoop
%   In this example, the Pantograph equation is not accretive and the method diverges.
%   However, by setting the option 'accretive=false', the Pantograph toolbox
%   solves the the anti-symmetrized system instead, which is accretive.
%
close all;

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(10/opt.pixel_size);
opt.boundaries_width = 800;
opt.termination_condition = TerminationCondition(tolerance= 1E-10, divergence_limit = 1E100);
opt.termination_condition_interval = 1;
opt.V_max = 0.5;
%opt.callback = DisplayCallback(plot_residual = true);
opt.forward_operator = true;
opt.crop_to_roi = false;

%% Medium parameters
lambda = 0.9;
a = 0.1 * ones(opt.N,1);
b = -5 * ones(opt.N,1);
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation. 
% First run the standard version, which will diverge because the
% equation is not accretive.
opt.accretive = true; % pretend that the equation is accretive, even though it is not
sim = PantographF(a, b, lambda, t0, opt);

% Define source
z = sim.grid.coordinates(1);
f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2);% - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z(1:t0));
%src = zeros(t0, 1);
%src(end) = 1;

source = sim.define_source(src);

%% Perform the different simulations and compare the results
precond = compare_simulations(sim, source, default_simulations);
no_precond = compare_simulations(sim, source, default_simulations, preconditioned=false);

%% Repeat with (anti-)symmetrized equation
s_opt = opt;
s_opt.accretive = false; % symmetrize system
s_sim = PantographF(a, b, lambda, t0, s_opt);
s_source = s_sim.define_source(src);

s_precond = compare_simulations(s_sim, s_source, default_simulations);

%%
%signal = [src; precond(end).value];
%s_signal = [src; s_precond(end).value];
%full_z = [z(1:t0); z+z(t0+1)];
%plot(full_z, real(signal));
%hold on;
%plot(full_z, real(s_signal));
%xlabel('t [s]');
%ylabel('x(t)');
%legend('Anysim symmetrized', 'BiCGStab non-symmetrized');

%%
[signal, state] = sim.exec(source);
[s_signal, s_state] = s_sim.exec(s_source);

%%
figure;
range = 1:200;
semilogy(state.residual_its(range), state.residuals(range), LineWidth = 2);
hold on;
semilogy(s_state.residual_its, s_state.residuals, LineWidth = 2);
hold off;
set(gca, FontSize=16);
xlabel('iteration');
ylabel('residual');
xlim([1, 200]);
ylim([1E-9, 1E18]);
yticks(10.^(-6:6:18));
legend(["original", "anti-symmetrized"], Location = "northwest", Box="off");
print(gcf, '-dpdf', 'converge_diverge_pantograph.pdf');
%%
%inspect_sim(sim);

