classdef Medium
    %MEDIUM Abstract base class for grid-based media
    % Used internally by simulation objects to prepare:
    % * the scattering potential V =  Tl (V_raw-V0) Tr, scaled to have 
    %   max ||V|| <= OPT.V_max. V_max should be below 1, and defaults
    %   to 0.95.
    %
    % * the (unscaled) background potential V0
    %
    % * the scaling matrices Tl and Tr used to scale the complete
    %   equation (see README.MD). In particular L = Tl (L_raw+V0) Tr
    %
    %   (c) 2021. Ivo Vellekoop
    properties
        G  % potential array = 1-V = 1-Tl*(Vraw-V0)*Tr
        V0 % background potential
        Tl % left scaling matrix
        Tr % right scaling matrix
        centers
        radii
        grid
        V_max % maximum norm for scaled V
    end
    
    methods
        function obj = Medium(V_raw, V_raw_min, grid, opt)
            % Construct the medium object
            %
            % This base class for TensorMedium and Diagonal medium
            % performs all common actions: input argument validation
            % and computation of V0.
            %
            % For Tensor media, V_raw should be reshaped so that all matrix
            % elements are in a column.
            %
            % The 'grid' object is used to add (tapered) boundaries are
            % to the scattering potential array.
            % 
            % Note: when the variations in V_raw are small, the scaling
            % matrices need to be very large to have max ||V|| == OPT.V_max
            % the values of V_raw - V0. If the matrices are too large, however
            % in the Green's function (L+1)^-1, the contribution of the +1
            % term will be relatively small. As a consequence, the Green's
            % function will have a large spatial extent, which may be larger
            % than the size of the simulation domain, causing wrapping artefacts.
            % To prevent this problem, V_raw_min is specified. Effectiely,
            % the algorithm pretends that the variations in V_raw-V0 are
            % always at least as large as V_raw_min. See analyze_dimensions
            % on how V_raw_min is calculated.
            %
            
            %% Set default options and validate input arguments
            defaults.V_max = 0.95;
            opt = set_defaults(defaults, opt);
            validateattributes(opt.V_max, {'numeric'}, {'scalar', '>', 0, '<', 1}); 
            obj.V_max = opt.V_max;

            sz = size(V_raw, 2+(1:grid.N_dim));
            if any(sz ~= grid.N_roi & sz ~= 1)
                error('Incorrect size for potential array');
            end
            
            %% Prepare scaling data by solving the smallest circle problem for 
            %% each individual component.
            % Note that this is optimal for scalar potentials and diagonal 
            % tensor potentials, but not necessarily optimal for full
            % tensor fields.
            % Note that at the moment only real elements are supported!!!
            N = size(V_raw, 1);
            M = size(V_raw, 2);
            centers = zeros(N, M);
            radii = zeros(N, M);
            for n=1:N
                for m=1:M
                    [centers(n,m), radii(n,m)] = smallest_circle(V_raw(n, m, :));
                end
            end
            obj.centers = centers;
            
            % If some elements of the potential matrix are (near) constant, 
            % some radii will be very small or zero, which affects precision
            % (and may cause wrap-around artefacts in some simulations). 
            % Therefore, enforce a minimum value that is used as radius
            %
            if any(radii < V_raw_min)
                warning('Slowing down simulation to prevent wrap-around artefacts')
            end
            obj.radii = max(radii, V_raw_min);
            obj.grid = grid;
        end
        
        function u = mix_source(obj, u, source, state)
            % MEDIUM.MIX_SOURCE(U, SOURCE, STATE) implements the
            % function G U + SOURCE, with G = 1-V and V the
            % scattering potential. 
            % SOURCE is a Source object (see Source)
            % STATE is a State object (see State), that is not
            % used in this implementation.
            u = source.apply(obj.multiplyG(u), state);
        end
        function u = mix_field(obj, u, uprop, state)
            % MEDIUM.MIX_FIELD(U, UPROP, STATE) implements the
            % function U + G (UPROP - U), with G = 1-V and V the
            % STATE is a State object, that is notified
            % about updates to the solution (see State).

            %u = u + obj.G .* (uprop - u);
            uprop = obj.multiplyG(uprop - u);
            if state.needs_report % checks termination condition
                M = norm(reshape(obj.grid.crop(uprop, 2), 1, []))^2;
                state.report_diff(M);
            end
            u = u + uprop;
        end 
    end
    methods (Abstract, Access=protected)
        u = multiplyG(obj, u)
    end
end

