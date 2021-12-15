classdef DilationMedium < ScalarMedium
    %DILATIONMEDIUM Helper class for the pantograph equation
    %   Implements (Δα + β Λs) f
    %   (c) 2021. Ivo Vellekoop
    properties
        coords  % s(t) t
        beta
    end
    methods
        function obj = DilationMedium(alpha, beta, s, grid, opt)
            % center coefficient 'α', put the constant part in 'L'
            obj@ScalarMedium(alpha, grid, opt);
            obj.coords = 1+(0:length(V)-1).* s;
            obj.beta = beta;
        end   
        function u = multiplyG(obj, u)
            u = obj.G .* u + obj.beta .* interp1(u, obj.coords);
        end  
    end
end

