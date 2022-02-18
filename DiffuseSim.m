classdef DiffuseSim < GridSim
    %DIFFUSESIM Solves the diffusion equation for inhomogeneous media
    %   Built on the AnySim framework.
    %   Solves either dynamic or steady state diffusion equation.
    %   Position-dependent absorption
    %   Position-dependent diffusion tensor or scalar coefficient
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = DiffuseSim(D, a, opt)
            % DIFFUSESIM Simulation object for a solving the diffusion
            % equation.
            %
            % sim = DiffuseSim(D, MUA, OPT) contructs a new simulation object
            % with the specified diffusion coefficient D and absoption
            % coefficient MUA.
            %
            % Note: in the sample code, for simplicity D and MUA are
            % specified in [um] and [um^-1], respectively, and micrometers
            % are used as unit for the temporal dimension.
            % Conversion to the regular units is done by scaling with  
            % the wave velocity: D -> D c and t -> t/c.
            %
            % Note: we always use 4 components (Fx, Fy, Fz, I). For
            % simulations in 1 or 2 dimensions, this approach is 
            % slightly wasteful since Fz and/or Fx will be 0.
            %
            % D Diffusion coefficient or 3x3 diffusion tensor.
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
            %                   [5 'um', 10 'um', 5 'um']. (default 1 '-')
            %
            % todo: allow indexed description for D (use an index to look up 
            %   the D tensor for each voxel).
            
            %% Set defaults
            opt.real_signal = true;
            defaults.pixel_size = 1;
            defaults.pixel_unit = 'm';
            defaults.V_max = 0.75; 
            defaults.alpha = 1;% theoretical optimum for Hermitian operators
            defaults.potential_type = "scalar"; 
            
            %defaults.V_max = 0.7208; %(sqrt(10)-1)/3
            %defaults.alpha = 2/(1+0.7208);
            
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            obj = obj@GridSim(4, opt); 
            
            %% Construct components: operators for medium, propagator and transform
            obj.makeMedium(D, a);
            obj.makePropagator();
        end
    end
    methods (Access = protected)        
        function makeMedium(obj, D, a)
            % Construct medium operator G=1-V
            %
            %    [      0]    
            %V = [  Q   0]
            %    [      0]
            %    [0 0 0 a]
            %
            %
            validatecompatiblesize(size(a), obj.opt.N);
            
            % Compute minimim values for max Re V
            % these values should be reachable to achieve sufficient
            % damping at the absorbing boundaries
            % mu = sqrt(Vraw(4) * Vraw)
            % how to choose Vraw(4)  ? 
            Vmin(4) = max(obj.mu_min);
            Vmin(1:3) = obj.mu_min(1:3).^2 / Vmin(4);

            % Invert D, then combines D and a into a single matrix
            % If the matrices are diagonal, store them as columns
            % instead of full matrices.
            D = data_array(D, obj.opt); % convert D to data array (put on gpu if needed, change precision if needed)
            szD = size(D);
            if obj.opt.potential_type == "tensor"
                % combine scalar 'a' and 3x3 matrix 'D' to 4x4 matrix
                validateattributes(D, {'numeric'}, {'nrows', 3, 'ncols', 3});
                validatecompatiblesize(szD(3:end), obj.opt.N);
                V = pageinv(D);
                V(4,4,:,:,:,:) = 0; 
                V = obj.grid.pad(V + padarray(shiftdim(a, -2), [3,3,0,0,0,0], 'pre'), 2);
                Vmin = diag(Vmin);
            elseif obj.opt.potential_type == "diagonal"
                % combine scalar 'a' and diagonal matrix 'D' to 4-element
                % diagonal matrix
                validateattributes(D, {'numeric'}, {'nrows', 3});
                validatecompatiblesize(szD(2:end), obj.opt.N);
                V = 1./D;
                V(4,:,:,:,:) = 0;
                V = obj.grid.pad(V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre'), 1);
            elseif obj.opt.potential_type == "scalar" 
                % combine scalar 'a' and scalar 'D' to 4-element diagonal
                % matrix
                obj.opt.potential_type = "diagonal";
                validatecompatiblesize(szD, obj.opt.N);
                V = repmat(shiftdim(1./D, -1), 3, 1, 1, 1, 1);
                V(4,:,:,:,:) = 0;
                V = obj.grid.pad(V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre'), 1);
            else
                error('Incorrect option for potential_type');
            end

            [obj.Tl, obj.Tr, obj.V0, V] = center_scale(V, Vmin, obj.opt.V_max);
            
            % apply scaling
            B = data_array(1 - V, obj.opt);
            if (obj.opt.potential_type == "diagonal")
                obj.medium = @(x) B .* x;
                obj.V0 = diag(obj.V0); % convert V0 to full matrix
            else
                obj.medium = @(x) fieldmultiply(B, x);
            end
        end
        
        function makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L the diffusion equation
            % differential operator.
            %
            % [               dx]
            % [     Q0        dy]
            % [               dz]
            % [dx dy dz (dt + a0*c)]
            V0 = obj.V0; % scaled background potential
            Tl = obj.Tl; % scaling matrices for the
            Tr = obj.Tr; % L operator
            
            Lr = zero_array([4, 4, obj.grid.N], obj.opt);
            
            % construct matrix with ikx, iky, ikz and iω
            for d=1:4
                location = zeros(4,4);
                location(4, d) = 1.0i;
                location(d, 4) = 1.0i;
                dd = location .* shiftdim(obj.grid.coordinates_f(d), -2);
                Lr = Lr + dd;
            end
            
            % L' + 1 = Tl (L+V0) Tr + 1
            % simplified to: L' = Tl L Tr + Tl V0 Tr + 1
            Lr = pagemtimes(pagemtimes(Tl, Lr), Tr) + (Tl * V0 * Tr + eye(4));
                        
            % invert L to obtain dampened Green's operator
            % then make x,y,z,t dimensions hermitian (set kx,ky,kz,omega to 0)
            % to avoid artefacts when N is even.
            % Note that there is a difference between first
            % taking the inverse and then taking the real partthen setting hermiIf we still need the forward operator Lr, 
            % we have to take special care at the the edges.
            % 
            Lr = pageinv(Lr);
            Lr = SimGrid.fix_edges_hermitian(Lr, 3:6);
            if obj.opt.forward_operator
                % note: this is very inefficient, we only need to do this
                % at the edges. However, the forward operator is only
                % needed for debugging purposes and for showing that
                % without preconditioner other methods (bicgstab, etc.) 
                % perform badly. So we don't care about optimizing this.
                LL = pageinv(Lr)-eye(4);
                obj.L = @(x) real(ifftv(fieldmultiply(LL, fftv(x))));
            end
            
            % the propagator just performs a
            % page-wise matrix-vector multiplication in k-space
            obj.propagator = @(x) real(ifftv(fieldmultiply(Lr, fftv(x))));
        end
    end
end
