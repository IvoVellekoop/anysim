classdef DilationMedium < DiagonalMedium
    %DILATIONMEDIUM Helper class for the pantograph equation
    %   Implements (Δα + β Λs) f
    %   (c) 2021. Ivo Vellekoop
    properties
        coords  % s(t) t
        beta
        s
    end
    methods
        function obj = DilationMedium(alpha, beta, s, grid, opt)
            % center coefficient 'α', put the constant part in 'L'
            if (isscalar(beta))
                beta = ones(1,1,grid.N(1)) * beta; 
            else
                beta = shiftdim(beta(:), -2);
            end
            if (isscalar(alpha))
                alpha = ones(1,1,grid.N(1)) * alpha; 
            else
                alpha = shiftdim(alpha(:), -2);
            end
            beta(1:opt.dilation_start) = 0;
            obj@DiagonalMedium(alpha, grid, opt);
            obj.coords = 1+(0:grid.N(1)-1).* s;
            obj.beta = beta * obj.Tl * obj.Tr;
            obj.s = s;
        end   
        function u = multiplyG(obj, u)
            scaled = reshape(interp1(u(:), obj.coords, 'linear', 0), size(u));
            u = obj.G .* u - obj.beta .* scaled;
        end  
    end
end

