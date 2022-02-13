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
            defaults.wavelength = 1;
            defaults.pixel_size = 0.25;
            defaults.pixel_unit = 'λ';
            defaults.N = size(n);
            defaults.V_max = 0.95;
            defaults.alpha = 0.75;%real(1/(1 + 1.0i*defaults.V_max));
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            % todo: we now always use 4 components (Fx, Fy, Fz, I)
            % for 2-dimensional simulations, Fz = 0, so this approach
            % is slightly wasteful.
            obj = obj@GridSim([], opt); 
            obj.k0 = 2*pi/opt.wavelength;

            %% Construct components: operators for medium, propagator and transform
            % Compute scaling factors. Note: Vraw = -i ε k0² = -i n² k0²            
            V = -1i * obj.k0^2 * n.^2;

            % Towards the edges, B is windowed, so that B=0 and V=1.
            % Since V = Tr .* (Vraw - V0) .* Tl, B0 implies that Vraw = V0 + 1/(Tr * Tl)
            % where V0 ~ -i k0^2
            % At this point, Vraw = V0 + Compute the minimum real part of that should be V.
            Vmin = imag((obj.k0 + 1i * max(obj.mu_min))^2); % minimum required absorption coefficient
            [obj.Tl, obj.Tr, obj.V0, V] = center_scale(V, Vmin, obj.opt.V_max);

            % apply scaling
            B = obj.grid.pad(data_array(1 - V, obj.opt), 0);
            obj.medium = @(x) B .* x;
            obj.Tl = obj.Tl * 1i; % include factor i to rotate source term??
            
            obj.propagator = obj.makePropagator();
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L+1)^-1 = 
            % with L = Tl(Lraw + V0)Tr, and Lraw the 
            % differential operator for the Helmholtz equation.
            % In k-space, Lraw = i (kx^2+ky^2+kz^2)
            
            % Compute -∇²=‖p‖² in k-space   
            % L' + 1 = [Tl (L+V0) Tr + 1]^-1
            L = obj.grid.coordinates_f(1).^2 + obj.grid.coordinates_f(2).^2 + obj.grid.coordinates_f(3).^2;
            L = (obj.Tl * obj.Tr) * (L - 1i * obj.V0); % -1i compensates for factor i in Tl matrix
            if obj.opt.forward_operator
                obj.L = @(u) ifftn(L .* fftn(u));
            end
            Lr = 1./(1+L);
            
            % point-wise multiplication in the Fourier domain
            propagator = @(u, state) ifftn(Lr .* fftn(u));
        end
     end
end
