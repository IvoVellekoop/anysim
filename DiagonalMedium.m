classdef DiagonalMedium
    %DIAGONALMEDIUM Helper class to implement a field potential of
    %diagonal tensors
    %
    %   (c) 2019. Ivo Vellekoop
    properties
        G  % potential array = 1-V = 1-Tl*(Vraw-V0)*Tr
        V0 % background potential
        Tl % left scaling matrix
        Tr % right scaling matrix
    end
    methods
        function obj = DiagonalMedium(V_raw, V_raw_min, grid, opt)
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
            opt = set_defaults(struct('V_max', 0.95), opt);
            validateattributes(opt.V_max, {'numeric'}, {'scalar', '>', 0, '<', 1}); 
            sz = size(V_raw, 1+(1:grid.N_dim));
            if any(sz ~= grid.N_roi & sz ~= 1)
                error('Incorrect size for potential array');
            end
            
            %% 1. Apply centering and scaling
            % For real diagonal matrices, this step is simple
            Nc = size(V_raw, 1);
            centers = zeros(Nc, 1);
            radii = zeros(Nc, 1);
            for r=1:Nc
                [centers(r), radii(r)] = smallest_circle(V_raw(r, :));
            end
            
            %% 2. Equilibrate matrix so that all elements have magnitude <= 1
            % If some elements of the potential matrix are (near) constant, 
            % some radii will be very small or zero, which affects precision
            % (and may cause wrap-around artefacts in some simulations). 
            % Therefore, enforce a minimum value that is used as radius
            %
            if any(radii < V_raw_min)
                warning('Slowing down simulation to prevent wrap-around artefacts')
            end
            radii = max(radii, V_raw_min);
            
            %% 3. Scale to have operator norm 1
            % Apply centering and scaling
            scale = opt.V_max./radii;
            obj.V0 = diag(centers);
            obj.Tl = diag(sqrt(scale));
            obj.Tr = diag(sqrt(scale));
            V = scale .* (V_raw - centers); 
            obj.G = grid.pad(data_array(1- V, opt), 1);
        end   
          
        function u = mix_source(obj, u, source, state)
            % MEDIUM.MIX_SOURCE(U, SOURCE, STATE) implements the
            % function G U + SOURCE, with G = 1-V and V the
            % scattering potential. 
            % SOURCE is a Source object (see Source)
            % STATE is a State object (see State), that is not
            % used in this implementation.
            u = source.apply(obj.G .* u, state);
        end
        function u = mix_field(obj, u, uprop, state)
            % MEDIUM.MIX_FIELD(U, UPROP, STATE) implements the
            % function U + G (UPROP - U), with G = 1-V and V the
            % STATE is a State object (see State), that is not
            % used in this implementation.

            %u = u + obj.G .* (uprop - u);
            uprop = obj.G .* (uprop - u);
            state.store_diff(uprop);
            u = u + uprop;
        end   
    end
end

