%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%
% This test simulates steady state diffusion in a 2-D medium with a ring
% with an anisotropic diffusion tensor (fast diffusion in tangential direction,
% slow in radial direction)
%

%% Set up simulation options
opt = DiffuseSimOptions();
opt.N = [256, 256];       % number of grid points in x,y
opt.boundaries_width = 0; % all boundaries periodic
opt.pixel_size = 0.25;
opt.pixel_unit = 'mm';
opt.callback = DisplayCallback(component = 4);
opt.forward_operator = true;

%% D is an anisotropic diffusion coefficient [10 0; 0 1], rotated over some angle
x = shiftdim(((1:opt.N(1))-opt.N(1)/2) / opt.N(1), -1);
y = shiftdim(((1:opt.N(2))-opt.N(2)/2) / opt.N(2), -2);
r = sqrt(x.^2 + y.^2);
x = x ./ r;
y = y ./ r;
tx = x;
x = y;
y = tx;

% define anisotropic ring diffusion coefficient
D = zeros([3,3,opt.N]);
D(1,1,:,:) = 10 * x.^2 + y.^2;
D(2,2,:,:) = x.^2 + 10 * y.^2;
D(3,3,:,:) = 1;
D(1,2,:,:) = -9 * x .* y;
D(2,1,:,:) = D(1,2,:,:);

% outside/inside ring: just low D
mask = r<0.2 | r > 0.3;
D(1,1, mask) = 0.5;
D(2,2, mask) = 0.5;
D(1,2, mask) = 0;
D(2,1, mask) = 0;

% absorbing boundary at the right (prevents wrap-around)
a = zeros(opt.N);
a(:,end-15:end) = 1;
sim = DiffuseSim(D, a, opt);

%% Place a source at the left
s = zeros(opt.N(1), opt.N(2));
s(:,1) = 1;
source = sim.define_source(shiftdim(s, -1), 4);
[u, state] = sim.exec(source);

%% Plot data
figure;
imagesc(squeeze(u(4,:,:)).');
hold on;
[Y,X] = meshgrid(1:opt.N(1), 1:opt.N(2));
U = squeeze(u(2,:,:));
V = squeeze(u(1,:,:));
startX = X(1,10:10:end);
startY = Y(1,10:10:end);
%streamline(Y, X, V, U, startY, startX);
r = 1:10:256;
quiver(Y(r,r), X(r,r), V(r,r), U(r,r), 5, 'k');

rectangle("Position",[0.2 * opt.N(2), 0.2* opt.N(1), 0.6*opt.N(2), 0.6*opt.N(1)], 'Curvature', 1);
rectangle("Position",[0.3 * opt.N(2), 0.3* opt.N(1), 0.4*opt.N(2), 0.4*opt.N(1)], 'Curvature', 1);
axis image;
xlim([16 240]);
ylim([16 240]);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
hold off;

%% Compare accuracies between simulation methods
sims = default_simulations("symmetric", has_adjoint = true);
%bare = compare_simulations(sim, source, default_simulations, preconditioned = false);
results = compare_simulations(sim, source, sims);
