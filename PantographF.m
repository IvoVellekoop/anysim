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
            opt.N_components = 2;
            obj = obj@GridSim(opt.N, opt); 

            %% Construct components: operators for medium, propagator and transform
            obj.t0 = t0;
            obj.lambda = lambda;
            obj = obj.makeMedium(alpha, beta);
            obj = obj.makePropagator();            
        end

        function u = finalize(obj, u)
            u = finalize@GridSim(obj, u);
            u = u(1, :).';
            u(obj.t0-1) = u(obj.t0-1) * 2; %compensate for finite pixel effect
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
            S = zeros(2, obj.grid.N);
            S(2, :) = interp1(values, coordinates, 'linear', 0);
            
            % the boundary condition at t0 is converted to a delta source
            S(2, 1) = S(2, 1) + values(end) / obj.grid.pixel_size;
        end
    end
    methods (Access = protected)        
        function obj = makeMedium(obj, alpha, beta)
            % Construct medium operator B=1-V
            % V includes the non-constant part of alpha, as well as the
            % effect of beta

            ar = min(real(alpha(:)));
            if (max(abs(beta(:))) > ar)
                warning('accretivity of the system is not guaranteed');
            end

            % note: we cannot use center_scale, because we also have to
            % account for the effect of beta.
            [obj.V0, alpha_radius] = smallest_circle(alpha);
            beta_radius = max(abs(beta(:)));
            obj.Tr = min(obj.opt.V_max/(alpha_radius + beta_radius), 1E3 * ar);
            obj.Tl = 1;

            alpha = (obj.data_array(alpha(:)) - obj.V0) * obj.Tr;
            beta = obj.data_array(beta(:).') * obj.Tr * sqrt(obj.lambda); %includes scaling factor of Λ
            beta_adj = conj(beta) / obj.lambda;
            
            coordinates = ((0:obj.grid.N - 1) + obj.t0) .* obj.lambda - obj.t0 + 1;
            coordinates_adj = ((0:obj.grid.N - 1) + obj.t0) ./ obj.lambda - obj.t0 + 1;
            alpha = shiftdim(alpha, -2);
            B = [0 1; 0 0] .* conj(alpha) + [0 0; -1 0] .* alpha + [1 0; 0 1];
            obj.medium = @(u) fieldmultiply(B, u) + [...
                interp1(beta_adj .* u(2,:), coordinates_adj, 'linear', 0);...
                -beta .* interp1(u(1,:), coordinates, 'linear', 0)];
        end

        function obj = makePropagator(obj)
            % L = Tl (i p + V0) Tr
            L = obj.Tr * (1i * shiftdim(obj.grid.coordinates_f(1), -2) + obj.V0);
            
            % construct L+1 inverse matrix:
            L1 = ([0 1; 0 0] .* conj(L) + [0 0; -1 0] .* L + [1 0; 0 1]) ./ (1+abs(L).^2);
            
            %Lh = obj.grid.fix_edges_hermitian(L1, 3);
            Lh = L1;
%            Lh(:,:,end/2) = 0;
            obj.propagator = @(u) ifftv(fieldmultiply(Lh, fftv(u)));
        end
    end
end
