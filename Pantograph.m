classdef Pantograph < GridSim
    %PANTOGRAPH Solves the pantograph equation with inhomogeneous coefficients
    %   Built on the AnySim framework.
    %   (-∂t + α + β Λ)x = b    for t >= t0
    %   x = b                   for t < t0
    %
    %   (c) 2021. Ivo Vellekoop & Tom Vettenburg
    %
    properties
        TrInternal
    end
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
            obj = obj@GridSim(opt.N, opt); 

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
            obj.Tr = 1;
            obj.TrInternal = min(obj.opt.V_max/(alpha_radius + beta_radius), 1E8 * abs(ar));
            coordinates = 1+(0:obj.grid.N(1)-1).'.* lambda;
            if (isscalar(alpha))
                alpha = alpha * ones(obj.grid.N(1), 1);
            end
            if (isscalar(beta))
                beta = beta * ones(obj.grid.N(1), 1);
            end

            alpha = (obj.data_array(alpha) - obj.V0) * obj.TrInternal;
            beta = obj.data_array(beta) * obj.TrInternal * sqrt(lambda); %includes scaling factor of Λ
            beta(1:t0-1) = 0;  % the part < t0 is included in L, so 
            alpha(1:t0-1) = 0; % V = 0 (meaning alpha=beta=0)
            
            Ba = 1-alpha;
            obj.medium = @(u) Ba .* u - beta .* interp1(u, coordinates, 'linear', 0);
            
            coordinates_adj = 1+(0:obj.grid.N(1)-1).' ./ lambda;
            beta_adj = conj(beta) / lambda;
            obj.medium_adj = @(u) conj(Ba) .* u - beta_adj .* interp1(u, coordinates_adj, 'linear', 0);
        end

        function obj = makePropagator(obj, start)
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
            Tr = obj.TrInternal; % scaling factor
            dt = obj.grid.pixel_size;

            rate = -exp(-obj.V0 * dt) * (dt - 2 * Tr) / (dt + 2 * Tr);
            step = -4 * dt * Tr / (dt^2 - 4 * Tr^2);
            r = dt/(dt + 2 * Tr);
            disp([step / dt * Tr, -(log(rate)+dt/Tr)/obj.V0/dt, r/step]);
            S = Tr * (-exp(dt * obj.V0/2) + 4 * dt / (dt + 2 * Tr))/(dt - 2 * Tr);

            %S = S - step;
            obj.propagator = @(u) Pantograph.convolve(start, rate, step, u, rate*(S+step), step-r); 
            obj.propagator_adj = @(u) Pantograph.convolve_adj(start, conj(rate), step, u, step-r); 
        end
    end
    methods (Static)
        function u = convolve(start, rate, step, u, S, r)
            %%
            % Operator (L+1)^{-1}
            % Performs the iteration (per element):
            uorg = u(start:end);

            u(start) = S * u(start-1) + u(start) * step;
            for n=start+1:length(u)
                u(n) = rate * u(n-1) + u(n) * step;
            end
            u(start:end) = u(start:end) - uorg * r;
            u(1:start-1) = 0.5 * u(1:start-1);
        end
        function u = convolve_adj(start, rate_conj, step, u, r)
            %%
            % Operator (L^*+1)^{-1}
            % Performs the iteration (per element):
            uorg = u(start:end);
            u(end) = u(end) * rate_conj;
            for n=length(u)-1:-1:start
                u(n) = rate_conj * u(n+1) + u(n) * step;
            end
            u(start:end) = u(start:end) - uorg * r;
            u(1:start-1) = 0.5 * u(1:start-1);
        end
    end
end
