%   Simple test of the DiffuseSim toolbox
%   (c) 2021. Ivo Vellekoop
%
% This test simulates steady state diffusion in a 1-D homogeneous medium
% with absorption, and compares the result to the analytical solution
%

%% Set up simulation options
opt = struct();                 % clear any previous options
opt.N = [256, 256, 1, 1];         % number of grid points in x,y,z,t 
opt.boundaries.periodic = true; % all boundaries periodic
opt.pixel_size = 0.5;
opt.pixel_unit = 'um';
opt.callback.handle = @DisplayCallback;
opt.callback.cross_section = {4};
opt.potential_type = 'tensor';

%% Construct medium 
a = 0 * ones(opt.N);    % absorption coefficient [um^-1]
D1 = [5 0 0; 0 5 0; 0 0 1];   % background
D2 = [25 0 0; 0 1 0; 0 0 1];  % anisotropic
D3 = [25 0 0; 0 25 0; 0 0 25];  % high
D4 = [1 0 0; 0 1 0; 0 0 1];  % low
D5 = [0 10 0; -10 0 0; 0 0 1]; % 'chiral'

M = 64;
range = 1:M;
D = repmat(D1, [1,1,opt.N]);
D(:,:, range, range + 100) = repmat(D2, [1,1,M,M]);
D(:,:, range + M, range + 100) = repmat(D3, [1,1,M,M]);
D(:,:, range + 2*M, range + 100) = repmat(D4, [1,1,M,M]);
D(:,:, range + 3*M, range + 100) = repmat(D5, [1,1,M,M]);
a(:,end-10:end) = 1;

sim = DiffuseSim(D, a, opt);

%% Place a source in at the left, and sink to the right
s = zeros(opt.N(1), opt.N(2));
s(:,1) = 1;
source = sim.define_source(shiftdim(s,-1), 4); % todo: check in 'define_source' if shiftdim was performed correctly
[u, state] = sim.exec(source);

%%
ucrop = u(:, :, 132 + (-70:70));
I = squeeze(ucrop(4,:,:));
X = 4:8:size(ucrop, 2);
Y = 4:8:size(ucrop, 3);
udown = ucrop(:,X,Y);
figure;
imagesc(I);
hold on;
quiver(Y, X, squeeze(udown(2,:,:)), squeeze(udown(1,:,:)));
hold off;
axis image;