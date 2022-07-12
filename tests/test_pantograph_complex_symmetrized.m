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
b = b / sqrt(0.5);
t0 = round(1/opt.pixel_size); % first second is starting condition

%% Set up AnySim simulation
sim = PantographF(a, b, lambda, t0, opt);
z = sim.grid.coordinates(1);
zdil = z(t0:end);
zdil = zdil - zdil(1) + opt.pixel_size;
% Define source
f_init = @(t)  exp(-50 * (t(:)-0.5).^2);% - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
src = f_init(z);
src = src(1:t0);
source = sim.define_source(src(:));

%% execute for b < 0 and repeat for b > 0
xn = [src; sim.exec(source)];
simp = PantographF(a, -b, lambda, t0, opt);
sourcep = simp.define_source(src(:));
xp = [src; simp.exec(sourcep)];
xn = xn(1:length(z));
xp = xp(1:length(z));

%% plot results
figure;
value_range = [min([min(real(xn)), min(imag(xn)), min(real(xp)), min(imag(xp))]), ...
    max([max(real(xn)), max(imag(xn)), max(real(xp)), max(imag(xp))])];
plots(1) = fill([z(1), z(1), 1.0, 1.0], 1.10 * [value_range, value_range([end, 1])], [0.9, 0.9, 1], 'EdgeColor', 'none'); hold on;
plots(2) = plot(z, real(xn), Color = [0, 0.75, 0], LineWidth = 3);
plots(3) = plot(z, imag(xn), Color = [0, 0.75, 0], LineStyle = ':', LineWidth = 3);
plots(4) = plot(z, real(xp), Color = [1.0, 0, 0], LineWidth = 1);
plots(5) = plot(z, imag(xp), Color = [1.0, 0, 0], LineStyle = ':', LineWidth = 1);
set(gca, FontSize = 12, XTick = 0:ceil(max(z)), LineWidth = 2, FontSize = 12);

xlabel('t [s]', FontSize = 14);
ylabel('x', FontSize = 14);
legend(plots([2, 4]), 'b = -5 s^{-1}', 'b = +5 s^{-1}', FontSize = 14)
ylim(1.10 * value_range);
print(gcf, '-depsc', 'pantograph.eps');

%% compare different iterative methods (for b negative)
simulations = default_simulations();
[precond, table] = compare_simulations(sim, source, simulations);
