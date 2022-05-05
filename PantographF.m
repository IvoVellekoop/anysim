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
            % equation, fft-based.
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
            opt.N_components = 2;
            obj = obj@GridSim(opt.N, opt); 

            %% Construct components: operators for medium, propagator and transform
            obj.t0 = t0;
            obj = obj.makeMedium(alpha, beta, lambda);
            obj = obj.makePropagator();            
        end

        function u = finalize(obj, u)
            u = finalize@GridSim(obj, u);
            u = u(1, :).';
            u(obj.t0-1) = u(obj.t0-1) * 2; %compensate for finite pixel effect
        end

        function S = define_source(obj, values, position)
            arguments
                obj
                values
                position = []
            end
            % the source is placed in the second component of the vector, the
            % adjoint part has source zero.
            S = define_source@GridSim(obj, shiftdim(values, -1), [2 position]);

            % the boundary condition at t0 is converted to a delta source
            S(2, obj.t0-1) = S(2, obj.t0-1) / obj.grid.pixel_size;
        end
    end
    methods (Access = protected)        
        function obj = makeMedium(obj, alpha, beta, lambda)
            % Construct medium operator B=1-V
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
            obj.Tr = min(obj.opt.V_max/(alpha_radius + beta_radius), 1E3 * ar);
            obj.Tl = 1;

            alpha = (obj.data_array(alpha(:)) - obj.V0) * obj.Tr;
            beta = obj.data_array(beta(:)) * obj.Tr * sqrt(lambda); %includes scaling factor of Λ
            if (isscalar(alpha))
                alpha = alpha * ones(obj.grid.N(1), 1);
            end
            if (isscalar(beta))
                beta = beta * ones(obj.grid.N(1), 1);
            end

            beta(1:obj.t0-1) = 0;  % the part < t0 is included in L, so 
            alpha(1:obj.t0-1) = 0; % V = 0 (meaning alpha=beta=0)
            beta = beta.';
            beta_adj = conj(beta) / lambda;
            
            coordinates = 1+(0:obj.grid.N(1)-1).* lambda;
            coordinates_adj = 1+(0:obj.grid.N(1)-1) ./ lambda;
            alpha = shiftdim(alpha, -2);
            B = [0 1; 0 0] .* conj(alpha) + [0 0; -1 0] .* alpha + [1 0; 0 1];
            obj.medium = @(u) fieldmultiply(B, u) + [...
                beta_adj .* interp1(u(2,:), coordinates_adj, 'linear', 0);...
                -beta .* interp1(u(1,:), coordinates, 'linear', 0)];
        end

        function obj = makePropagator(obj)
            % L = Tl (i p + V0) Tr
            L = obj.Tr * (1i * shiftdim(obj.grid.coordinates_f(1), -2) + obj.V0);
            % construct L+1 inverse matrix:
            L1 = ([0 1; 0 0] .* conj(L) + [0 0; -1 0] .* L + [1 0; 0 1]) ./ (1+abs(L).^2);
            Lh = obj.grid.fix_edges_hermitian(L1, 3);
            start_matrix = inv([1, -obj.Tr; obj.Tr, 1]);
            obj.propagator = @(u) PantographF.convolve(obj.t0, Lh, u, start_matrix);
        end
    end
    methods (Static)
        function u = convolve(start, kernel, u, start_matrix)
            s = u(:, 1:start-2);
            u(:, 1:start-2) = 0;
            u = ifftv(fieldmultiply(kernel, fftv(u)));
            u(:, 1:start-2) = fieldmultiply(start_matrix, s);
        end
    end
end
