%   Simple test of the HelmholtzSim toolbox
%   (c) 2021. Ivo Vellekoop
%

%% Set up simulation options
opt = HelmholtzSimOptions(); % clear any previous options
opt.N = 256;
opt.boundaries_width = 64;
opt.callback = DisplayCallback();
opt.callback.interval = 100;
opt.termination_condition = TerminationCondition(relative_limit = 1E-6, iteration_count = 1E6);
opt.forward_operator = true; % for testing and comparison with MATLAB algorithms
opt.crop_to_roi = false;
opt.alpha = 0.7;%real(1/(1 + 1.0i*defaults.V_max));
lambda_r = 0.8:0.05:1;%1.2;


n = ones(opt.N); % refractive index
n(1:80) = repmat([0.9; 1.1], [40, 1]);
n(end-79:end) = repmat([0.9; 1.1], [40, 1]);

energies = zeros(length(lambda_r), 1);
iterations = zeros(length(lambda_r), 1);

for lambda_i = 1:length(lambda_r)
    opt.wavelength = lambda_r(lambda_i);
    sim = HelmholtzSim(n, opt);
                
    %% Define source and run the simulation
    source = sim.define_source(ones(opt.N(1)/2,1,1)); % intensity-only source (isotropic) at t=0
    [E, state] = sim.exec(source);
    energies(lambda_i) = norm(E(:))^2;
    iterations(lambda_i) = state.iteration;
end

%%
plot(energies);

[L, GL] = simulation_eigenvalues(sim);



