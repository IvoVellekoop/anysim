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
c = x ./ r; % cosine of angle
s = y ./ r; % sine of angle

% define anisotropic ring diffusion coefficient, rotated so that the fast
% axis of the diffusion tensor is in the tangential direction of the ring.
% 
% [c  s] [1  0 ] [c -s] = [c  s] [c     -s]  = [c^2 + 10s^2  9cs]
% [-s c] [0  10] [s  c] = [-s c] [10s  10c]  = [9cs          s^2 + 10c^2]


% [c  s] [1  0 ] [c -s] = [c  s] [10c -10s]  = [10c^2 + s^2  -9cs]
% [-s c] [0  10] [s  c] = [-s c] [  s    c]  = [-9cs         10 s^2 + c^2]
%
% for x=0 (top and bottom of ring, need fast diffusion along x direction)
% [10 0]
% [0  1]
D = zeros([3,3,opt.N]);
Dmax = 25;
Dmin = 1;
D(1,1,:,:) = Dmin * c.^2 + Dmax * s.^2;
D(2,2,:,:) = Dmax * c.^2 + Dmin * s.^2;
D(3,3,:,:) = Dmin;
D(1,2,:,:) = (Dmax-Dmin) * c .* s;
D(2,1,:,:) = D(1,2,:,:);

% outside/inside ring: just low D
mask = r<0.2 | r > 0.3;
D(1,1, mask) = 2;
D(2,2, mask) = 2;
D(3,3, mask) = 2;
D(1,2, mask) = 0;
D(2,1, mask) = 0;

% absorbing boundary at the bottom (prevents wrap-around)
a = zeros(opt.N);
a(:, end-15:end) = 1;
sim = DiffuseSim(D, a, opt);

%% Place a source at the top
s = zeros(opt.N(1), opt.N(2));
s(:, 1) = 1;
%s(128,end-70) = 5;
source = sim.define_source(shiftdim(s, -1), 4);
[u, state] = sim.exec(source);

%% Plot data
figure;
x = sim.grid.coordinates(1);
y = sim.grid.coordinates(2);
imagesc(x, y, squeeze(u(4,:,:)).');
hold on;
[X, Y] = meshgrid(x, y);
U = squeeze(u(1,:,:));
V = squeeze(u(2,:,:));
r = 16:16:240;
%r = 16:5:240;
quiver(X(r,r), Y(r,r), U(r,r).', V(r,r).', 1, 'k');

rectangle("Position",[0.2 * x(end), 0.2* y(end), 0.6 * x(end), 0.6 * y(end)], 'Curvature', 1);
rectangle("Position",[0.3 * x(end), 0.3* y(end), 0.4 * x(end), 0.4 * y(end)], 'Curvature', 1);
axis image;
xlim([x(16) x(240)]);
ylim([y(16) y(240)]);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
hold off;
print(gcf, '-depsc', 'diffusion_ring.eps');


%% Compare accuracies between simulation methods
sims = default_simulations(has_adjoint = true);
%bare = compare_simulations(sim, source, default_simulations, preconditioned = false);
[results, table] = compare_simulations(sim, source, sims);
