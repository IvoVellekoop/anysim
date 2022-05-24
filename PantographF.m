classdef PantographF < GridSim
    %PANTOGRAPH Solves the pantograph equation with inhomogeneous coefficients
    %   Built on the AnySim framework.
    %   (-∂t + α + β Λ)x = b    for t >= t0
    %   x = b                   for t < t0
    %
    %   (c) 2021. Ivo Vellekoop & Tom Vettenburg
    %
    properties
        t0
        lambda
        beta
    end
    methods
        function obj = PantographF(alpha, beta, lambda, t0, opt)
            arguments
                alpha (:,1)
                beta (:,1)
                lambda (1,1) {mustBePositive}
                t0 (1,1) {mustBeInteger, mustBePositive}
                opt PantographOptions
            end
            % PANTOGRAPHF Simulation object for a solving the pantograph
            % equation, fft-based
            % 
            % sim = PANTOGRAPHF(ALPHA, BETA, LAMBDA, T0, OPT)
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
            if opt.non_accretive
                opt.N_components = 2;
            else
                opt.N_components = 1;
            end
            obj = obj@GridSim(opt.N, opt); 

            %% Construct components: operators for medium, propagator and transform
            obj.t0 = t0;
            obj.lambda = lambda;
            obj = obj.makeMedium(alpha, beta);
            obj = obj.makePropagator();            
        end

        function u = finalize(obj, u)
            u = finalize@GridSim(obj, u).';
            if obj.opt.non_accretive
                u = u(:, 1);
            end
        end

        function S = define_source(obj, values)
            arguments
                obj
                values (:, 1)
            end
            % Constructs the source term Λ x_0
            % adds a delta at the first element to take into account the
            % boundary condition.
            if length(values) ~= obj.t0
                error("Source data must have size %d", obj.t0);
            end

            % rescale by lambda. Note: simulation starts at t0
            coordinates = (obj.t0 + (0:obj.grid.N-1)) * obj.lambda;
            S = obj.beta(:).' .* interp1(values, coordinates(:).', 'linear', 0) * (obj.Tl/obj.Tr);
            
            % the boundary condition at t0 is converted to a delta source
            S(1, 1) = S(1, 1) + values(end) / obj.grid.pixel_size;

            if obj.opt.non_accretive
                S = [zeros(1, obj.grid.N); S];
            end
        end
    end
    methods (Access = protected)        
        function obj = makeMedium(obj, alpha, beta)
            % Construct medium operator B=1-V
            % V includes the non-constant part of alpha, as well as the
            % effect of beta

            % note: we cannot use center_scale, because we also have to
            % account for the effect of beta.
            [obj.V0, alpha_radius] = smallest_circle(alpha);
            beta_radius = max(abs(beta(:)));
            obj.Tr = min(obj.opt.V_max/(alpha_radius + beta_radius), 1E3 * abs(min(real(alpha))));
            obj.Tl = 1;

            alpha = (obj.data_array(alpha(:)) - obj.V0) * obj.Tr;
            obj.beta = obj.data_array(beta(:).') * obj.Tr; %includes scaling factor of Λ
            
            % construct a sparse 'sampling' matrix for scaling the signal
            before = 1:obj.grid.N;
            after = (before - 1 + obj.t0) .* obj.lambda - obj.t0 + 1;
            fa = floor(after);
            ca = ceil(after);
            N = obj.grid.N;
            mask_f = (fa <= N) & (fa > 1);
            mask_c = (ca <= N) & (ca > 1);
            s1 = sparse(before(mask_f), fa(mask_f), fa(mask_f) - after(mask_f) + 1, N, N);
            s2 = sparse(before(mask_c), ca(mask_c), 1 + after(mask_c) - ca(mask_c), N, N);
            sample = (obj.beta(:).' .* max(s1, s2)).';
            
            if obj.opt.non_accretive
                alpha = shiftdim(alpha, -2);
                B = [0 1; 0 0] .* conj(alpha) + [0 0; -1 0] .* alpha + [1 0; 0 1];
                obj.medium = @(u) fieldmultiply(B, u) + [u(2,:) * sample'; -u(1,:) * sample];
            else
                B = 1 - alpha.';
                obj.medium = @(u) B.*u - u * sample;
            end
        end

        function obj = makePropagator(obj)
            % constructs the Green's function
            % there are two ways to do this: 
            rate = obj.V0 + 1/obj.Tr;
            G = exp(-rate * obj.grid.coordinates(1)) * obj.grid.pixel_size/obj.Tr;
            %L = shiftdim(1./fft(G) - 1, -2);
            
            L = obj.Tr * (1i * obj.grid.coordinates_f(1) + obj.V0);
            Linv = 1./(L+1);
            Lginv = fft(G);
            Lcomb = real(Linv) + 1i * imag(Lginv);
            L = 1./Lcomb-1;
            
            
            % construct L+1 inverse matrix:
            %filt = shiftdim(fftshift(tukeywin(obj.grid.N, 0.1)), -2);
            if obj.opt.non_accretive
                L = shiftdim(L, -2);
                Lh = ([0 1; 0 0] .* conj(L) + [0 0; -1 0] .* L + [1 0; 0 1]) ./ (1+abs(L).^2);
                obj.propagator = @(u) ifftv(fieldmultiply(Lh, fftv(u)));
            else
                obj.propagator = @(u) ifft(Lcomb.' .* fft(u));
            end
        end
    end
end
