classdef Pantograph < GridSim
    %PANTOGRAPH Solves the pantograph equation with inhomogeneous coefficients
    %   Built on the AnySim framework.
    %   (- ∂t f + α + β Λs f + εΘ(t0−t))
    %
    %   (c) 2021. Ivo Vellekoop & Tom Vettenburg
    %
    properties
        r_dil   % norm of V_raw, excluding the variations in a
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
            defaults.pixel_size = {1, 's'};
            defaults.alpha = 0.5;
            defaults.V_max = 0.999; % theoretical optimum for Hermitian operators
            opt = set_defaults(defaults, opt);
            
            %% Construct base class
            obj = obj@GridSim([], opt); 
            
            %% Construct components: operators for medium, propagator and transform
            % note: change minus sign for alpha and beta, so that positive
            % alpha corresponds to absorption (consistent with A)
            obj.r_dil = max(abs(beta(:))) / sqrt(s);
            obj.medium  = obj.makeMedium(-alpha, -beta, s);
            obj.transform  = NoTransform();
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
            % L + 1 = Tl (dt + V0) Tr + 1
            % Now construct (L+1)^-1 u = 
            % (conv(u, exp(-α t)) + u(t0) exp(-α t)) / (Tl Tr)  
            % with α = V0 + 1/Tl Tr
            % approximate finite difference: each dx = pixel-size
            % step, decrease signal by exp(-α dx)
            V0 = obj.medium.V0; % scaled background potential
            s = obj.medium.Tl * obj.medium.Tr; % scaling factor
            alpha = V0 + 1 / s;
            rate = exp(-alpha * obj.grid.pixel_size(1));
            start = obj.opt.dilation_start;
            propagator.apply = @(u, state) Pantograph.convolve(start, rate, obj.grid.pixel_size(1), s/(s+1), u) / s; 
        end
        
        function [centers, radii, feature_size, bclimited] = adjustScale(obj, centers, radii)
            % 
            % The current radii are only based on variations in 'a'.
            % We need to also account for the dillation part
            radii = radii + obj.r_dil;
            bclimited = false;
            feature_size = 1/abs(centers+radii);
        end
    end
    methods (Static)
        function u = convolve(start, rate, step, f, u)
            u(1:start-1) = u(1:start-1) * f;
            for n=start:length(u)
                u(n) = rate * u(n-1) + u(n) * step;
            end
        end
    end
end
