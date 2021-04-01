classdef DiagonalMedium < Medium
    %DIAGONALMEDIUM Helper class to implement a field potential of
    %diagonal tensors
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = DiagonalMedium(V_raw, V_raw_min, grid, opt)
            obj@Medium(V_raw, V_raw_min, grid, opt);
            
            %% Scale to have operator norm OPT.V_max
            scale = obj.V_max./obj.radii;
            obj.V0 = diag(obj.centers);
            obj.Tl = eye(length(scale));
            obj.Tr = diag(scale);
            V = scale .* (V_raw - obj.centers); 
            obj.G = grid.pad(data_array(1- V, opt), 2);
        end   
    end
    methods (Access=protected)
        function u = multiplyG(obj, u)
            u = obj.G .* u;
        end  
    end
end

