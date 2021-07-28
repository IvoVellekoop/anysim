classdef HelmholtzSim < GridSim
    %HELMHOLTZSIM Solves the Helmholtz equation for inhomogeneous media
    %   Built on the AnySim framework.
    %
    %   (c) 2021. Ivo Vellekoop
    properties
        k0
    end
    methods
        function obj = HelmholtzSim(n, opt)
            % HELMHOLTZSIM Simulation object for a solving the Helmholtz
            % equation.
            %
            % After rotation, the equation is given by:
            % -i ∇² u + ⟨V⟩ u + ΔV u = i S 
            % with V = -i ε k0²
            % sim = DiffuseSim(N, K0, OPT) contructs a new simulation object
            % with the specified refractive index N.
            %
            % n: refractive index distribution. 
            % The positive imaginary part corresponds to absorption.
            %
            % Options:
            % OPT.pixel_size    pixel pitch, may be different in different
            %                   dimensions (not recommended).
            %                   Recommended: set actual SI units ( e.g.
            %                   0.25 'μm' )
            %                   (default 0.25 [-])
            % OPT.wavelength    wavelength. (default 1).
            % OPT.N             Number of pixels in the simulation.
            %                   (defaults to size(N), no need to set
            %                   explictly unless using singleton expansion)
            %
            %% Set defaults
            opt.real_signal = false;
            opt.potential_type = "scalar";
            defaults.wavelength = 1;
            defaults.pixel_size = 0.25;
            defaults.pixel_unit = 'λ';
            defaults.N = size(n);
            defaults.V_max = 0.618034; % optimal for V covering full unit disk. Todo: improve further?
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            % todo: we now always use 4 components (Fx, Fy, Fz, I)
            % for 2-dimensional simulations, Fz = 0, so this approach
            % is slightly wasteful.
            obj = obj@GridSim([], opt); 
            
            %% Construct components: operators for medium, propagator and transform
            % Note: V = -i ε k0² = -i n² k0²
            obj.k0 = 2*pi/opt.wavelength;
            obj.medium  = DiagonalMedium(-1i * obj.k0^2 * n.^2, obj.grid, obj.opt);
            obj.medium.Tl = obj.medium.Tl * 1i; % include factor i to rotate source term
            obj.transform  = FourierTransform(obj.opt);
            obj.propagator = obj.makePropagator();
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L the 
            % differential operator for the Helmholtz equation.
            % In k-space, we simply have L= i (kx^2+ky^2+kz^2)
            V0 = -1i * obj.medium.V0; % raw (non-scaled) background potential. Compensate for factor i in Tl matrix
            Tl = obj.medium.Tl; % scaling factor, contains factor i
            Tr = obj.medium.Tr; % scaling factor (real)
            
            % Compute -∇²=‖p‖² in k-space   
            % L' + 1 = [Tl (L+V0) Tr + 1]^-1
            Lr = obj.grid.coordinates_f(1).^2 + obj.grid.coordinates_f(2).^2 + obj.grid.coordinates_f(3).^2;
            Lr = obj.to_internal(Tl*Tr*(Lr + V0));
            if obj.opt.forward_operator
                obj.L = @(u) pagemtimes (Lr, u);
            end
            Lr = 1./(1+Lr);
            
            % point-wise multiplication
            propagator.apply = @(u, state) Lr .* u;
        end
    end
    methods (Access=protected)
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
