classdef HelmholtzSim < GridSim
    %HELMHOLTZSIM Solves the Helmholtz equation for inhomogeneous media
    %   Built on the AnySim framework.
    %
    %   (c) 2021. Ivo Vellekoop
    properties (SetAccess = private)
        k0
    end
    methods
        function obj = HelmholtzSim(n, opt)
            arguments
                n double
                opt HelmholtzSimOptions
            end
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
            
            %% Construct base class
            opt = opt.validate(size(n));
            obj = obj@GridSim(opt.N, opt); 
            obj.k0 = 2*pi/opt.wavelength;

            %% Construct components: operators for medium, propagator and transform
            obj = obj.makeMedium(n);
            obj = obj.makePropagator();
        end
    end
    methods (Access = protected)        
        function obj = makeMedium(obj, n)
            % Compute scaling factors. Note: Vraw = -i ε k0² = -i n² k0²            
            V = -1i * obj.k0^2 * n.^2;

            % Towards the edges, B is windowed, so that B=0 and V=1.
            % Since V = Tr .* (Vraw - V0) .* Tl, B0 implies that Vraw = V0 + 1/(Tr * Tl)
            % where V0 ~ -i k0^2
            % At this point, Vraw = V0 + Compute the minimum real part of that should be V.
            Vmin = imag((obj.k0 + 1i * max(obj.mu_min))^2); % minimum required absorption coefficient

            if obj.opt.legacy_mode
                [V1, V2] = bounds(imag(V(:)));
                obj.V0 = 1i * (V1 + V2) / 2;
                V = V - obj.V0;
                obj.Tl = 1;
                obj.Tr = obj.opt.V_max/max(abs(V(:)));
                V = obj.Tr * V;
            else
                [obj.Tl, obj.Tr, obj.V0, V] = center_scale(V, Vmin, obj.opt.V_max);
            end

            % apply scaling
            B = obj.grid.pad(obj.data_array(1 - V), 0);
            obj.medium = @(x) B .* x;
            obj.medium_adj = @(x) conj(B) .* x;
            obj.Tl = obj.Tl * 1i; % include factor i to rotate source term??
        end        
        function obj = makePropagator(obj)
            % Constructs the propagator (L+1)^-1 = 
            % with L = Tl(Lraw + V0)Tr, and Lraw the 
            % differential operator for the Helmholtz equation.
            % In k-space, Lraw = i (kx^2+ky^2+kz^2)
            
            % Compute -∇²=‖p‖² in k-space   
            % L' + 1 = [Tl (L+V0) Tr + 1]^-1
            L = obj.grid.coordinates_f(1).^2;
            for d = 2:obj.grid.N_dim
                L = L + obj.grid.coordinates_f(d).^2;
            end
            
            L = (obj.Tl * obj.Tr) * (L - 1i * obj.V0); % -1i compensates for factor i in Tl matrix
            if obj.opt.forward_operator
                obj.L = @(u) ifftn(L .* fftn(u));
            end
            Lr = 1./(1+L);
            
            % point-wise multiplication in the Fourier domain
            obj.propagator = @(u, state) ifftn(Lr .* fftn(u));
            obj.propagator_adj = @(u, state) ifftn(conj(Lr) .* fftn(u));
        end
     end
end
