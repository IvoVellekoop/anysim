classdef DiagonalMedium < Medium
    %DIAGONALMEDIUM Helper class to implement a field potential of
    %diagonal tensors
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = DiagonalMedium(V_raw, V_raw_min, grid, opt)
            obj@Medium(V_raw, V_raw_min, grid, opt);
            
            %% Scale to have operator norm OPT.V_max
            scale = opt.V_max./obj.radii;
            obj.V0 = diag(obj.centers);
            obj.Tl = eye(4);
            obj.Tr = diag(scale);
            V = scale .* (V_raw - obj.centers); 
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

