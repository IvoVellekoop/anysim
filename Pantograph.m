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
            obj.Tr = min(obj.opt.V_max/(alpha_radius + beta_radius), 1E8 * ar);
            coordinates = 1+(0:obj.grid.N(1)-1).'.* lambda;
            if (isscalar(alpha))
                alpha = alpha * ones(obj.grid.N(1), 1);
            end
            if (isscalar(beta))
                beta = beta * ones(obj.grid.N(1), 1);
            end

            alpha = (obj.data_array(alpha) - obj.V0) * obj.Tr;
            beta = obj.data_array(beta) * obj.Tr;
            beta(1:t0-1) = 0;  % the part < t0 is included in L, so 
            alpha(1:t0-1) = 0; % V = 0 (meaning alpha=beta=0)
            
            Ba = 1-alpha;
            obj.medium = @(u) Ba .* u - beta .* interp1(u, coordinates, 'linear', 0);
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
            Tr = obj.Tl * obj.Tr; % scaling factor
            dt = obj.grid.pixel_size;
            F1 = 1/Tr;%-log(1-dt/Tr) / dt; % approx 1/Tr
            rate = exp(-(obj.V0 + F1) * dt);
            step = exp(-(obj.V0 + 0.5*F1) * dt) * dt / Tr;
            disp([F1 * Tr, step / dt * Tr]);
 
            F = 1/(Tr+1); % multiplication factor for u0
            S = 0;%step/2;% subtraction factor for uort
            obj.propagator = @(u) Pantograph.convolve(start, rate, step, 1/(Tr+1), u, F, S); 
        end
    end
    methods (Static)
        function u = convolve(start, rate, step, f, u, F, S)
            %%
            % Operator (L+1)^{-1}
            % Performs the iteration (per element):
            % $u[n] \rightarrow c_4 u[n]$               for $n<start$
            % $u[n] \rightarrow c_3 u[n-1] + c_2 u[n]$  for $n=start$
            % $u[n] \rightarrow c_1 u[n-1] + c_2 u[n]$  for $n>start$
            % 
            % followed by multiplication by c0
            % 
            % Result for full vector (take start = 0 for simplicity)
            %
            % $u[n] \rightarrow c_04 u[n]$                           for $n<0$
            % $u[n] \rightarrow c_034 u[-1] c_1^n + 
            %       sum_{t=0}^n c_02 u[t] c_1^{n-t}$                 for $n>=0$
            %
            % c_04 = c_0 c_4, etc.
            % Solution is at:
            % B[1-(L+1)^{-1}B] u = B(L+1)^{-1} b 
            % [1-(L+1)^{-1}B] u = (L+1)^{-1} b 
            % u-(L+1)^{-1}(B u + b) = 0
            %
            % Take B scalar const for n>=0 and 1 for n<0 so that solution is exponential decay.
            %
            % source: (b[n] = 0 for n >= 0)
            % (L+1)^{-1} b = c_04 b[n]                      for $n<0$
            %              = c_034 b[-1] c_1^n              for $n>=0$
            %
            % solution: 
            % u[n] = b[n]           for $n<0$
            %      = b[-1] q^(n+1)  for $n>=0$
            %
            % (L+1)^{-1} u[n] = 
            %      = c_04 b[n]            for $n<0$
            %      = c_034 b[-1] c_1^n + 
            %        q sum_{t=0}^n c_02 b[-1] q^t c_1^{n-t}$                 for $n>=0$
            %
            % simplify for n>=0
            %      = b[-1] c_1^n (c_034  + c_02 q sum_{t=0}^n (q/c_1)^t)
            %      = b[-1] c_1^n (c_034  + c_02 q (1-(q/c_1)^{n+1})/(1-q/c_1)
            %      = c_034 b[-1] c_1^n   + c_02 q b[-1] (c_1^{n+1}-q^{n+1})/(c_1-q)
            %
            % substitute:  u-(L+1)^{-1}(B u + b) = 0
            % for n < 0
            % b[n] - 2 c_04 b[n] = 0          => c_04 = 1/2
            %
            % for n >=0
            % b[-1] q^(n+1)
            %  - B (c_034 b[-1] c_1^n + c_02 q b[-1] (c_1^{n+1}-q^{n+1})/(c_1-q))
            %  - c_034 b[-1] c_1^n = 0
            % 
            % B = 0 -> q = c_1
            % b[-1] q^(n+1) - c_034 b[-1] c_1^n = 0
            % c_034 = c_1 
            %
            % Required: c_1 <= q  (L+1)^{-1} more damped than (L+V)^{-1}
            % equate c_1^n terms for nonzero B:
            %  - B (c_034 + c_02 q c_1/(c_1-q)) - c_034 = 0
            %  - B (1 + c_02 q /(c_1-q)) - 1 = 0
            %  - B ((c_1-q + c_02 q) /(c_1-q)) = 1
            %  (c_1-q + c_02 q) /(q-c_1) = 1/B
            % 
            % equate q^(n+1) terms:
            %  1 + B (c_02 q /(c_1-q)) = 0
            %  c_02 q /(q-c1) = 1/B
            %
            % combine:
            % (c_1-q)/(q-c_1) = 0!?
            u0 = u(start-1);
            u(start-1) = u(start-1) * F;% * 1.46;%e1.333333;
            uorg = u(start:end);
            for n=start:length(u)
                u(n) = rate * u(n-1) + u(n) * step;
            end
            u(start:end) = u(start:end) - uorg * S;
            u(start-1) = u0;
            u(1:start-1) = u(1:start-1) * f;
        end
    end
end
