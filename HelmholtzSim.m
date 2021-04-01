classdef HelmholtzSim < GridSim
    %HELMHOLTZSIM Solves the Helmholtz equation for inhomogeneous media
    %   Built on the AnySim framework.
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = HelmholtzSim(n, opt)
            % HELMHOLTZSIM Simulation object for a solving the diffusion
            % equation.
            %
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
            obj.medium  = obj.makeMedium(1i * (n * 2*pi / opt.wavelength).^2);
            obj.medium.Tl = obj.medium.Tl * 1i; % rotate source. Todo: sign??
            obj.transform  = FourierTransform(obj.opt);
            obj.propagator = obj.makePropagator();
        end
        
        function propagator = makePropagator(obj)
            % Constructs the propagator (L'+1)^-1 = 
            % with L' = Tl(L + V0)Tr, and L the 
            % differential operator for the Helmholtz equation.
            % In k-space, we simply have L=-(kx^2+ky^2+kz^2)
            V0 = obj.medium.V0; % scaled background potential
            Tl = obj.medium.Tl; % scaling factors for the
            Tr = obj.medium.Tr; % L operator
            
            Lr = obj.grid.coordinates_f(1).^2 + obj.grid.coordinates_f(2).^2 + obj.grid.coordinates_f(3).^2;
            
            % L' + 1 = [Tl (L+V0) Tr + 1]^-1
%            Lr = obj.to_internal((1 + Tl*Tr*V0 - 1i*Tl*Tr*Lr).^(-1));
            Lr = obj.to_internal((1 + Tl*Tr*(Lr + 1.0i * V0)).^(-1));
            
            % point-wise multiplication
            propagator.apply = @(u, state) Lr .* u;
        end
    end
    methods (Access=protected)
        function Vmin = analyzeDimensions(obj, Vmax)
            
            %% The Green's function [L+1]^-1 decays exponentially
            %% with a decay coefficient of 
            Vmin = 1
        end
    end
end
