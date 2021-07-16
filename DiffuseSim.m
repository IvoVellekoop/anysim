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
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            obj = obj@GridSim(4, opt); 
            
            %% Construct components: operators for medium, propagator and transform
            obj.medium  = obj.makeMedium(D, a);
            obj.transform  = FourierTransform(obj.opt);
            obj.propagator = obj.makePropagator();
        end
    end
    methods (Access = protected)        
        function medium = makeMedium(obj, D, a)
            % Construct medium operator G=1-V
            %
            %    [      0]    
            %V = [  Q   0]
            %    [      0]
            %    [0 0 0 a]
            %
            %% Invert D, then add entries for a
            D = data_array(D, obj.opt); % convert D to data array (put on gpu if needed, change precision if needed)
            if obj.opt.potential_type == "tensor"
                % combine scalar 'a' and 3x3 matrix 'D' to 4x4 matrix
                validateattributes(D, {'numeric'}, {'nrows', 3, 'ncols', 3});
                V = pagefun(@inv, D);
                V(4,4,:,:,:,:) = 0; 
                V = V + padarray(shiftdim(a, -2), [3,3,0,0,0,0], 'pre');
                medium = TensorMedium(obj.to_internal(V), obj.grid, obj.opt);
            elseif obj.opt.potential_type == "diagonal"
                % combine scalar 'a' and diagonal matrix 'D' to 4-element
                % diagonal matrix
                validateattributes(D, {'numeric'}, {'nrows', 3});
                V = 1./D;
                V(4,:,:,:,:) = 0;
                V = V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre');
                medium = DiagonalMedium(obj.to_internal(V), obj.grid, obj.opt);
            elseif obj.opt.potential_type == "scalar" 
                % combine scalar 'a' and scalar 'D' to 4-element diagonal
                % matrix
                obj.opt.potential_type = "diagonal";
                V = repmat(shiftdim(1./D, -1), 3, 1, 1, 1, 1);
                V(4,:,:,:,:) = 0;
                V = V + padarray(shiftdim(a, -1), [3,0,0,0,0], 'pre');
                medium = DiagonalMedium(obj.to_internal(V), obj.grid, obj.opt);
            else
                error('Incorrect option for potential_type');
            end
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L the diffusion equation
            % differential operator.
            %
            % [               dx]
            % [     Q0        dy]
            % [               dz]
            % [dx dy dz (dt + a0*c)]
            V0 = obj.medium.V0; % scaled background potential
            Tl = obj.medium.Tl; % scaling matrices for the
            Tr = obj.medium.Tr; % L operator
            
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
            Lr = pagemtimes(pagemtimes(Tl, Lr + V0), Tr);
                        
            % invert L to obtain dampened Green's operator
            % then make x,y,z,t dimensions hermitian (set kx,ky,kz,omega to 0)
            % to avoid artefacts when N is even.
            % Note that there is a difference between first
            % taking the inverse and then taking the real partthen setting hermiIf we still need the forward operator Lr, 
            % we have to take special care at the the edges.
            % 
            Lr = pageinv(Lr + eye(4));
            Lr = SimGrid.fix_edges_hermitian(Lr, 3:6);
            if obj.opt.forward_operator
                % note: this is very inefficient, we only need to do this
                % at the edges. However, the forward operator is only
                % needed for debugging purposes and for showing that
                % without preconditioner other methods (bicgstab, etc.) 
                % perform badly. So we don't care about optimizing this.
                LL = pageinv(Lr)-eye(4);
                obj.L = @(u) pagemtimes (LL, u);
            end
            
            % the propagator just performs a
            % page-wise matrix-vector multiplication
            propagator.apply = @(u, state) pagemtimes(Lr, u);
        end
        
        function [centers, radii, feature_size, bclimited] = adjustScale(obj, centers, radii)
            % For the diffusion equation, the attenuation coefficient is
            % given by sqrt(mu_a Q), so we have a choice whether to
            % increase mu_a or Q. We want our choice to have as least 
            % impact on the convergence rate as possible.
            % To 'guess' an optimum, we choose the combination for which
            % the decay along each dimensions is proportional
            % to the size of the simulation domain.
            
            % determine required absorption coefficients
            % so that solution decays sufficiently fast
            if isvector(centers)
                mu_min = obj.mu_min;
                V_raw_max = centers + radii;
                mu_current = sqrt(V_raw_max(4) * V_raw_max);% current maximum absorption coefficients

                if any(mu_current < mu_min)
                    % todo: not necessarily the optimal splitting between Vmin(4)
                    %       and Vmin(1:3)
                    % todo: untested code! test for mu_min > 0
                    
                    V_max_2(4) = max(mu_min(4), norm(V_max(1:3)));
                    V_max_2(1:3) = mu_min(1:3).^2 / V_max_2(4);
                    shift = V_max_2 - V_raw_max;
                    centers = centers + shift/2;
                    radii = radii + shift/2;
                    bclimited = true;
                    warning('untested code!');
                else
                    % no need to adjust
                    bclimited = false;                    
                end
                feature_size = 1./sqrt(V_raw_max * V_raw_max(4)); %todo: correct for steady state?
            else
                % for tensors: only process diagonal
                [centers, radii, feature_size, bclimited] = diag(obj.scale_adjuster(diag(centers), diag(radii)));
            end
        end
    end
end
