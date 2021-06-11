classdef ScalarMedium < Medium
    %SCALARMEDIUM Helper class to implement a scalar scattering potential
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = ScalarMedium(V_raw, V_raw_min, grid, opt)
            obj@Medium(V_raw, V_raw_min, grid, opt);
            
            %% Scale to have operator norm OPT.V_max
            scale = obj.V_max/obj.radii;
            obj.Tl = 1;
            obj.Tr = scale;
            obj.V0 = obj.centers;
            V = scale * (V_raw - obj.centers); 
            obj.G = grid.pad(data_array(1 - V, opt), 2);
        end   
        function u = multiplyG(obj, u)
            u = obj.G .* u;
        end  
    end
end

