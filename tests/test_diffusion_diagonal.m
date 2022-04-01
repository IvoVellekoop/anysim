%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%
% This test simulates steady state diffusion in a 1-D homogeneous medium
% with absorption, and compares the result to the analytical solution
%

%% Set up simulation options
opt = DiffuseSimOptions();                 % clear any previous options
opt.N = [256, 256];
opt.boundaries_width = 0; % all boundaries periodic
opt.pixel_size = 0.5;
opt.pixel_unit = 'um';
opt.callback = DisplayCallback(component = 4);
opt.forward_operator = true;

%% Construct medium 
a = 0.01;    % absorption coefficient [um^-1]
D1 = [5; 5; 1];   % background
D2 = [25; 1; 1];  % anisotropic
D3 = [25; 25; 25];  % high
D4 = [1; 1; 1];  % low
D5 = D2;%[0 10 0; 10 0 0; 0 0 1]; % 'chiral'

M = 64;
range = 1:M;
D = repmat(D1, [1,opt.N]);
D(:, range, range + 100) = repmat(D2, [1,M,M]);
D(:, range + M, range + 100) = repmat(D3, [1,M,M]);
D(:, range + 2*M, range + 100) = repmat(D4, [1,M,M]);
D(:, range + 3*M, range + 100) = repmat(D5, [1,M,M]);

sim = DiffuseSim(D, a, opt);

%% Place a source in at the left
s = zeros(opt.N(1), opt.N(2));
s(:,1) = 1;
source = sim.define_source(shiftdim(s,-1), 4); % todo: check in 'define_source' if shiftdim was performed correctly
[u, state] = sim.exec(source);
I = squeeze(u(4,:, :));

%%
X = 5:10:opt.N(1);
Y = 5:10:opt.N(2);
udown = u(:,X,Y);
figure;
imagesc(I);
hold on;
quiver(X, Y, squeeze(udown(2,:,:)), squeeze(udown(1,:,:)));
hold off;
axis image;

%% Compare to other methods and compute errors
simulations = default_simulations(has_adjoint = true);

% without preconditioner, all methods diverge!
% bare = compare_simulations(sim, source, simulations, preconditioned = false);
[precond, table] = compare_simulations(sim, source, simulations);
