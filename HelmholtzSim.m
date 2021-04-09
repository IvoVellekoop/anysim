classdef HelmholtzSim < GridSim
    %HELMHOLTZSIM Solves the Helmholtz equation for inhomogeneous media
    %   Built on the AnySim framework.
    %
    %   (c) 2021. Ivo Vellekoop
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
            defaults.pixel_size = {0.25, '-'};
            defaults.N = size(n);
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            % todo: we now always use 4 components (Fx, Fy, Fz, I)
            % for 2-dimensional simulations, Fz = 0, so this approach
            % is slightly wasteful.
            obj = obj@GridSim([], opt); 
            
            %% Construct components: operators for medium, propagator and transform
            % Note: V = -i ε k0² = -i n² k0²
            obj.medium  = obj.makeMedium(-1i * (n * 2*pi / opt.wavelength).^2);
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
            Lr = obj.grid.coordinates_f(1).^2 + obj.grid.coordinates_f(2).^2 + obj.grid.coordinates_f(3).^2;
            
            % L' + 1 = [Tl (L+V0) Tr + 1]^-1
            Lr = obj.to_internal((1 + Tl*Tr*(Lr + V0)).^(-1));
            
            % point-wise multiplication
            propagator.apply = @(u, state) Lr .* u;
        end
    end
    methods (Access=protected)
        function Vmin = analyzeDimensions(obj, Vmax)
            
            % The Green's function [L+1]^-1 is given by
            % 1/(i T (‖p‖² + V0) + 1)
            % with T = Tr Tl / i the scaling coefficient and V0=-⟨ε⟩k0²
            % Performing an 1-D inverse Fourier transform, we find
            % that this function corresponds to a bidirectional decaying
            % exponential with exponent α=sqrt(V0 + i/T))
            % For low absorption, we can expand
            % α ≈ sqrt(V0) + i/(2 T sqrt(V0)) 
            %   = sqrt(V0) + 1/(2 T sqrt(-V0))
            % 
            % The damping is given by the real part of α
            % α = Re(sqrt(V0) + 1/(2 T sqrt(-V0)))
            %
            % So, the minimum value for T to ensure a decay that is fast
            % enough is given by
            % T >= 1/2 Re(1/ sqrt(-V0)) / (α-Re(sqrt(V0)))
            %
            limiting_size = obj.grid.dimensions();
            limiting_size = max(limiting_size(~obj.grid.boundaries.periodic));
            alpha_min = 10 ./ limiting_size;
            %Tmin = 
            sV1 = real(sqrt(Vmax/1i));    %Re(sqrt(V0)
            sV2 = real(1/sqrt(Vmax*1i));  %Re(1/ sqrt(-V0))
            Vmin = 0.5 * sV2 / (alpha_min-sV1);
        end
    end
end
