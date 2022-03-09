classdef Pantograph < GridSim
    %PANTOGRAPH Solves the pantograph equation with inhomogeneous coefficients
    %   Built on the AnySim framework.
    %   (-∂t + α + β Λ)x = b    for t >= t0
    %   x = b                   for t < t0
    %
    %   (c) 2021. Ivo Vellekoop & Tom Vettenburg
    %
    methods
        function obj = Pantograph(alpha, beta, lambda, t0, opt)
            arguments
                alpha (:,1)
                beta (:,1)
                lambda (1,1) {mustBePositive}
                t0 (1,1) {mustBeInteger, mustBePositive}
                opt PantographOptions
            end
            % PANTOGRAPH Simulation object for a solving the pantograph
            % equation
            % 
            % sim = PANTOGRAPH(ALPHA, BETA, LAMBDA, T0, OPT)
            % contructs a new simulation object with the specified
            % coefficients. 
            % T0 is the start time given in gridpoints.
            % Values of ALPHA(1:T0-1) and BETA(1:T0-1) are ignored.
            %
            % ∂t x = α x + β Λ x
            % with Λ the unitary dilation operator with scale factor λ:
            % Λ x(t) = sqrt(λ) x(λ t)
            
            %% Construct base class
            opt = opt.validate(size(alpha), size(beta));
            obj = obj@GridSim(opt.N, opt, opt); 

            %% Construct components: operators for medium, propagator and transform
            obj = obj.makeMedium(alpha, beta, t0, lambda);
            obj = obj.makePropagator(t0);            
        end
    end
    methods (Access = protected)        
        function obj = makeMedium(obj, alpha, beta, t0, lambda)
            % Construct medium operator G=1-V
            % V includes the non-constant part of alpha, as well as the
            % effect of beta

            ar = min(real(alpha(:)));
            if (max(abs(beta(:))) > ar)
                warning('acrretivity of the system is not guaranteed');
            end

            % note: we cannot use center_scale, because we also have to
            % account for the effect of beta.
            [obj.V0, alpha_radius] = smallest_circle(alpha);
            beta_radius = max(abs(beta(:)));
            obj.Tl = 1;
            obj.Tr = min(1/(alpha_radius + beta_radius), 1E8 * ar);

            coordinates = 1+(0:obj.grid.N(1)-1).'.* lambda;
            if (isscalar(alpha))
                alpha = alpha * ones(obj.grid.N(1), 1);
            end
            if (isscalar(beta))
                beta = beta * ones(obj.grid.N(1), 1);
            end
            
            alpha = (obj.data_array(alpha) - obj.V0) * obj.Tr;
            beta = obj.data_array(beta) * obj.Tr;
            beta(1:t0) = 0;  % the part < t0 is included in L+1, so 
            alpha(1:t0) = 0; % V = 0 (meaning alpha=beta=0)
            
            B = 1-alpha;
            obj.medium = @(u) B .* u(:)- beta .* interp1(u(:), coordinates, 'linear', 0);
        end

        function obj = makePropagator(obj, t0)
            % L + 1 = Tl (dt + V0) Tr + 1   for t >= t0, 
            % Tl Tr + 1                     for t < t0
            %
            % define s = Tl Tr
            % L + 1 = s [dt + (V0 + 1/s)]   for t >= t0
            % L + 1 = s + 1                 for t < t0
            %
            % define q = V0 + 1/s and construct (L+1)^-1 u
            % (conv(u, exp(-q t)) + u(t0) exp(-q t)) / s    for t >= t0
            % 1 / (s+1)                                     for t < t0
            %
            % implement as rolling convolution. Each grid step, decrease signal by exp(-q dt)
            % u'[t<t0] <= u[t<0] * s / (s + 1)
            % u'[t] = u'[t-1] * exp[-q dt]  +  dt u[t]
            %
            % In the last line, the last term converts each pixel in 'u'
            % to a 'source'. The last first term computes a rolling convolution.
            % 
            s = obj.Tl * obj.Tr; % scaling factor
            q = obj.V0 + 1 / s;
            rate = exp(-q * obj.grid.pixel_size(1));
 %           rate = 1 - q * obj.grid.pixel_size(1); % finite difference
            obj.propagator = @(u) Pantograph.convolve(t0, rate, obj.grid.pixel_size(1), s/(s+1), u) / s; 
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
