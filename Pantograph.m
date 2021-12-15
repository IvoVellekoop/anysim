classdef Pantograph < GridSim
    %PANTOGRAPH Solves the pantograph equation with inhomogeneous coefficients
    %   Built on the AnySim framework.
    %   (- ∂t f + α + β Λs f + εΘ(t0−t))
    %
    %   (c) 2021. Ivo Vellekoop & Tom Vettenburg
    %
    properties
        alpha, beta, s % coefficients of the differential equation
    end
    methods
        function obj = Pantograph(alpha, beta, s, opt)
            % PANTOGRAPH Simulation object for a solving the diffusion
            % equation.
            %
            % sim = Pantograph(alpha, beta, s, t0)
            % contructs a new simulation object with the specified
            % coefficients. t0 is the start time given in the same units
            % as pixel_size. Values of alpha, beta, s in the range 0-t0
            % will be ignored.
            %
            % OPT Options structure
            %   .pixel_size     Grid spacing, specified as, for example
            %                   [5 'ms']. (default 1 's')
            
            %% Set defaults
            opt.real_signal = true;
            defaults.pixel_size = {1, 's'};
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            obj = obj@GridSim(1, opt); 
            
            %% Construct components: operators for medium, propagator and transform
            obj.medium  = obj.makeMedium(alpha, beta, s);
            obj.transform  = FourierTransform(obj.opt);
            obj.propagator = obj.makePropagator();
        end
    end
    methods (Access = protected)        
        function medium = makeMedium(obj, alpha, beta, s)
            % Construct medium operator G=1-V
            % perform scaling so that ‖V‖ < 1
            medium = DilationMedium(alpha, beta, s, obj.grid, obj.opt);
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L = -iω  
            %
            V0 = obj.medium.V0; % scaled background potential
            Tl = obj.medium.Tl; % scaling matrices for the
            Tr = obj.medium.Tr; % L operator
            
            Lr = shiftdim(obj.grid.coordinates_f(1), -2);
            
            % L' + 1 = Tl (L+V0) Tr + 1
            % simplified to: L' = Tl L Tr + Tl V0 Tr + 1
            Lr = Tl * Tr * (Lr + V0) + 1;
                        
            % invert L to obtain dampened Green's operator
            % then make x,y,z,t dimensions hermitian (set kx,ky,kz,omega to 0)
            % to avoid artefacts when N is even.
            % Note that there is a difference between first
            % taking the inverse and then taking the real partthen setting hermiIf we still need the forward operator Lr, 
            % we have to take special care at the the edges.
            % 
            Lr = 1./ Lr;
            Lr = SimGrid.fix_edges_hermitian(Lr, 2);
            if obj.opt.forward_operator
                % note: this is very inefficient, we only need to do this
                % at the edges. However, the forward operator is only
                % needed for debugging purposes and for showing that
                % without preconditioner other methods (bicgstab, etc.) 
                % perform badly. So we don't care about optimizing this.
                LL = 1./Lr-1;
                obj.L = @(u) LL .* u;
            end
            
            % the propagator just performs a
            % point-wise multiplication
            propagator.apply = @(u, state) Lr .* u;
        end
        
        function [centers, radii, feature_size, bclimited] = adjustScale(obj, centers, radii)
            % 
            % We have V_raw = -i n^2 k0^2
            % n = sqrt(i V_raw) / k0
            %
            % The imaginary part of n corresponds to absorption
            % The boundaries tend towards G=0 -> V_raw=centers+radii,
            % so the absorption should be strong enough there
            V_max_abs = centers + radii;
            k1 = sqrt(1.0i * V_max_abs);
            mu_min = max(obj.mu_min); % minimum required absorption coefficient
            if imag(k1) < mu_min
                % replace absorption by minimum required absorption
                % and calculate what V_raw corresponds to this point
                k2 = real(k1) + 1.0i * mu_min;
                V_raw = -1.0i * k2^2; 
                
                % Adjust centers and radii so that new V_raw is included
                %todo: not exactly optimal!
                shift = V_max_abs - V_raw;
                centers = centers - shift/2;
                radii = radii + abs(shift)/2;
                bclimited = true;
            else
                bclimited = false;
            end
            
            %todo: not exact! Largest real part of n not necessarily reached at n0
            feature_size = pi / real(k1);            
        end
    end
end
