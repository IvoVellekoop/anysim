%   Simple test of the Pantograph toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Simulation parameters
opt = PantographOptions(); % clear any previous options
opt.pixel_size = 0.01;
opt.N = round(8/opt.pixel_size);
%opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
%opt.callback = DisplayCallback();
opt.boundaries_width = 0; % don't add boundaries

%% Medium parameters
lambda = 0.5;
a = zeros(opt.N, 1) + 5;
a(round(6/opt.pixel_size) : end) = 5 - 10i;  
b = zeros(opt.N, 1) - 5;
b(round(3/opt.pixel_size) : round(5)/opt.pixel_size) = 0;
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
sim = Pantograph(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
zdil = z(t0:end);
zdil = zdil - zdil(1) + opt.pixel_size;
% Define source
f_init = @(t)  exp(-50 * (t(:)-0.5).^2);% - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z);
src = src(1:t0-1);
source = sim.define_source(src(:));

%% execute for b < 0 and repeat for b > 0
xn = sim.exec(source);
simp = Pantograph(a, -b, lambda, t0, opt);
sourcep = simp.define_source(src(:));
xp = simp.exec(source);

%% plot results
figure;
hold on;
plot(z, real(xp), 'b', LineWidth = 2);
plot(z, real(xn), 'r', LineWidth = 1);
plot(z, imag(xp), 'b--', LineWidth = 2);
plot(z, imag(xn), 'r--', LineWidth = 1);
set(gca, FontSize = 12);
xlabel('t [s]', FontSize = 14);
ylabel('x', FontSize = 14);
legend('b = 5 s^{-1}', 'b = -5 s^{-1}', FontSize = 14)
print(gcf, '-depsc', 'pantograph.eps');

%% compare different iterative methods (for b negative)
simulations = default_simulations("nonsymmetric");
[precond, table] = compare_simulations(sim, source, simulations);
