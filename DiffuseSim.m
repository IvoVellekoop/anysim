classdef DiffuseSim < GridSim
    %DIFFUSESIM Solves the diffusion equation for inhomogeneous media
    %   Built on the AnySim framework.
    %   Solves either dynamic or steady state diffusion equation.
    %   Position-dependent absorption
    %   Position-dependent diffusion tensor or scalar coefficient
    %
    %   (c) 2019. Ivo Vellekoop
    methods
        function obj = DiffuseSim(D, a, opt)
            % DIFFUSESIM Simulation object for a solving the diffusion
            % equation.
            %
            % sim = DiffuseSim(D, MUA, OPT) contructs a new simulation object
            % with the specified diffusion coefficient D and absoption
            % coefficient MUA.
            %
            % Alternatively, D and MUA are be specified in [m] and [m^-1],
            % respectively, using meters as unit for the spatial dimension.
            % Conversion to the regular units is done by multiplying or 
            % dividing by the wave velocity: D -> D c and t -> t/c
            %
            % D Diffusion coefficient or diffusion tensor. 
            %   The dimensions of D depend on the OPT.potential_type setting.
            %   Assuming grid dimensions Nx x Ny x Nz x Nt the size must be:
            %   'scalar' size(D) = [Nx, Ny, Nz, Nt]
            %   'diagonal' size(D) = [3, Nx, Ny, Nz, Nt]
            %   'tensor' size (D) = [3, 3, Nx, Ny, Nz, Nt]
            %   In all cases, D must be strict positive definite.
            %   Singleton dimensions are expanded automatically.
            %
            % MUA Absorption coefficient.
            %   Must be a positive scalar (may be 0) of size [Nx, Ny, Nz,
            %   Nt]. Singleton dimensions are expanded automatically.
            %
            % OPT Options structure
            %   .potential_type Any of the options 'scalar' (default), 'diagonal', 'tensor'
            %                   This option determines how the D array is
            %                   interpreted (see above)
            %   .pixel_size     Grid spacing, specified as, for example
            %                   [5 'um', 10 'um', 5 'um']. Note: at the
            %                   moment, all dimensions should have the same unit!
            %
            %   todo: allow indexed description for D (use an index to look up 
            %   the D tensor for each voxel).
            %
            %   todo: allow specifying full 4 x 4 Qa tensor so that anisotropic
            %   absorption coefficients are possible.
            
            %% Set defaults
            opt.real_signal = true;
            
            %% Construct base class
            % todo: we now always use 4 components (Fx, Fy, Fz, I)
            % for 2-dimensional simulations, Fz = 0, so this approach
            % is slightly wasteful.
            obj = obj@GridSim(4, opt); 
            
            %% Construct components: operators for medium, propagator and transform
            obj.medium  = obj.makeMedium(D, a);
            obj.transform  = FourierTransform(obj.opt);
            obj.propagator = obj.makePropagator();
        end 
        
        function medium = makeMedium(obj, D, a)
            % Construct L and V operators
            %
            %    [      0]    
            %V = [  Q   0] => T ./ ([Dx, Dy, Dz, Dt]' * [Dx, Dy, Dz, Dt]) 
            %    [      0]
            %    [0 0 0 a]
            % For simplicity, we always convert the combination of a and Q to a 4x4 tensor field
            % todo: allow for a block-diagonal representation of the direct sum a x Q
            % todo: allow for a completely diagonal representation of a x Q 
            % todo: allow for a scalar representation of a x Q
            
            %% Invert D, then add entries for a
            D = data_array(D, obj.opt); % convert D to data array (put on gpu if needed, change precision if needed)
            if obj.opt.potential_type == "tensor"
                validateattributes(D, {'numeric'}, {'nrows', 3, 'ncols', 3});
                % invert 'D', then convert 'V' and 'a' to dimension
                % 4*4*(Nx or 1)*(Ny or 1)*(Nz or 1)*(Nt or 1)
                V = pagefun(@inv, D);
                V(4,4,:,:,:,:) = 0; 
                V = V + padarray(shiftdim(a, -2), [3,3,0,0,0,0], 'pre');                
            elseif obj.opt.potential_type == "diagonal"
                validateattributes(D, {'numeric'}, {'nrows', 3});
                % invert 'D', then convert 'V' and 'a' to dimension
                % 4*(Nx or 1)*(Ny or 1)*(Nz or 1)*(Nt or 1)
                V = 1./D;
                V(4,:,:,:,:) = 0;
                V = V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre');
            elseif obj.opt.potential_type == "scalar" %convert to diagonal matrix
                % first convert scalar 'D' to 3x3 diagonal, then
                % proceed as with diagonal 'D'
                obj.opt.potential_type = "diagonal";
                V = repmat(shiftdim(1./D, -1), 3, 1, 1, 1, 1);
                V(4,:,:,:,:) = 0;
                V = V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre');
            else
                error('Incorrect option for potential_type');
            end
            medium = makeMedium@GridSim(obj, V);
        end
        
        function Vmin = analyzeDimensions(obj, Vmax)
            % Analyzes the dimensions and spacings of the simulation
            % grid. Issues a warning if the grid spacing is larger
            % than the size of the smallest features.
            % Also computes a minimum value for the scattering potential
            % to ensure that the Green's function does not extend beyond
            % the size of the simulation window.
            
            % The scaled Green's function (L+1)^-1 decays exponentially in
            % space and time. This decay should be fast enough to ensure
            % convergence (this scaling is performed by makeMedium by
            % ensuring that the scaled potential V has ||V||<1).
            % The decay should also be fast enough to minimize wrap-around
            % artefacts. For this reason, we require the decay length/time
            % to be equal to a fraction of the size/duration of the
            % simulation in each dimension.
            % 
            % The decay coefficient for the diffusion equation is given by
            % mu_eff = sqrt(mu_a / D}        (assuming c=1)
            % Or, in term of elements of the un-scaled potential: 
            % mu_eff_j = sqrt(Vraw_t * Vraw_j)    
            % with j=x,y,z, or t. So, Vraw_t = mu_eff_t
            % and V_raw_j = mu_eff_j^2 / Vraw_t
            %
            % start by computing minimum required mu_eff in all dimensions

            active = obj.grid.N ~= 1;
            limiting_size = obj.grid.dimensions().';
            % requires same pixel size in all dimensions?
            limiting_size(obj.grid.boundaries.periodic) = max(limiting_size);
            mu_eff_min = 10 ./ limiting_size;
            
            % special case for steady state (t-axis inactive)
            if ~active(4)
                mu_eff_min(4) = max(mu_eff_min);
            end
            Vmin = mu_eff_min.^2 / mu_eff_min(4);
            
            % Now, check if the resolution is high enough to resolve the
            % smallest features
            Vmax = max(Vmin, Vmax);
            feature_size = 1./sqrt(Vmax(4) * Vmax(active)); %todo: correct for steady state?
            pixel_size = obj.grid.pixel_size(active);
            
            if any(feature_size/2 < pixel_size)
                res_limit = sprintf("%g ", feature_size/2);
                res_current = sprintf("%g ", pixel_size);
                warning("Resolution is too low to resolve the smallest features in the simulation. Minimum pixel size: [%s] Current pixel size: [%s]", res_limit, res_current);
            end
            if any(feature_size/8 > pixel_size)
                res_limit = sprintf("%g ", feature_size/2);
                res_current = sprintf("%g ", pixel_size);
                warning("Resolution seems to be on the high side. Minimum pixel size: [%s] Current pixel size: [%s]", res_limit, res_current);
            end
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L the diffusion equation
            % differential operator.
            % note: requires that grid and medium properties already set!
            % V0, Tl, Tr and grid are properties of medium 
            % [               dx]
            % [     Q0        dy]
            % [               dz]
            % [dx dy dz (dt + a0*c)]
            V0 = obj.medium.V0; % scaled background potential
            Tl = obj.medium.Tl; % scaling matrices for the
            Tr = obj.medium.Tr; % L operator
            
            Lr = zero_array([4, 4, obj.grid.N], obj.opt);
            
            % insert dx, dy, dz, dt
            for d=1:4
                location = zeros(4,4);
                location(4, d) = 1.0i;
                location(d, 4) = 1.0i;
                dd = location .* shiftdim(obj.grid.coordinates_f(d), -2);
                Lr = Lr + dd;
            end
            
            % L' + 1 = Tl (L+V0) Tr + 1
            % simplified to: L' = Tl L Tr + Tl V0 Tr + 1
            Lr = pagemtimes(pagemtimes(Tl, Lr), Tr);
            Lr = Lr + Tl * V0 * Tr + eye(4);            
            
            % invert L to obtain dampened Green's operator
            if isa(Lr, 'gpuArray') || isa(Lr, 'distributed')
                Lr = pagefun(@inv, Lr);
            else % why is pagefun not implemented for ordinary arrays?
                sLr = size(Lr);
                Lr = Lr(:,:,:);
                for n=1:size(Lr, 3)
                    Lr(:,:,n) = inv(Lr(:,:,n));
                end
                Lr = reshape(Lr, sLr);
            end
            
            Lr = SimGrid.fix_edges_hermitian(Lr, 3:6); % make x,y,z,t dimensions hermitian
            
            % page-wise matrix-vector multiplication
            propagator.apply = @(u, state) permute(...
                pagemtimes(Lr, permute(u, [1,6,2,3,4,5])),...
                [1,3,4,5,6,2]);
        end
    end
end
