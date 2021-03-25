classdef TensorMedium < Medium
    %TENSORMEDIUM Helper class to implement a tensor field potential
    %
    %   (c) 2021. Ivo Vellekoop
    methods
        function obj = TensorMedium(V_raw, V_raw_min, grid, opt)
            Nc = size(V_raw, 1);
            if size(V_raw, 2) ~= Nc
                error('Tensors in scattering potential must be square')
            end
            obj@Medium(reshape(V_raw, Nc * Nc, []), V_raw_min, grid, opt)
            
            %% Equilibrate matrix so that all elements have magnitude <= 1
            % note: this is not necessarily optimal, but it is a simple
            % way to (more or less) evenly scale the variance in all dimensions
            radii = reshape(obj.radii, Nc, Nc);
            centers = reshape(obj.centers, Nc, Nc);
            [P,R,C] = equilibrate(radii);
            obj.Tl = P'*R*P;
            obj.Tr = C;
            obj.V0 = centers;
            V = pagemtimes(pagemtimes(obj.Tl, V_raw - obj.V0), obj.Tr);
            
            %% The procedure above does not guarantee that ||V||< OPT.V_max
            % In this final step, we compute the norm of V everywhere
            % and adjust the global scaling as needed
            Vnorm = max(pagefun(@norm, V), [], 'all');
            scale = opt.V_max./Vnorm;
            obj.Tl = obj.Tl * sqrt(scale);
            obj.Tr = obj.Tr * sqrt(scale);
            
            obj.G = grid.pad(data_array(1- V, opt), 2);
        end
    end
    methods (Access=protected)
        function u = multiplyG(obj, u)
            u = pagemtimes(obj.G, u);
        end
    end
end

