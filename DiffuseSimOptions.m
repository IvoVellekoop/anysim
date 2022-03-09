classdef DiffuseSimOptions < AnySimOptions & GridOptions
    %DIFFUSESIMOPTIONS Options for simulation of the diffusion equation
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % Dimensions of the grid in voxels, excluding boundaries.
        % When empty, the size is determined automatically from D and a.
        N (1,:) {mustBePositive} = []
        potential_type string {mustBeMember(potential_type, ["scalar", "diagonal", "tensor"])} = []
    end
    methods 
        function obj = DiffuseSimOptions(opt)
            arguments 
                opt.?DiffuseSimOptions
            end
            obj.alpha = 1; % theoretical optimum for Hermitian operators
            obj = copy_properties(obj, opt);
        end
        function opt = validate(opt, szD, sza)
            % the diffusion simulation is always 3-dimensional, so we need
            % 3 flux components + 1 intensity
            opt.N_components = 4;

            % 'guess' potential type if it is not provided
            if isempty(opt.potential_type)
                if szD(1) == 3 && szD(2) == 3
                    opt.potential_type = "tensor";
                elseif szD(1) == 3
                    opt.potential_type = "diagonal";
                else                  
                    opt.potential_type = "scalar"; 
                end
            end
            
            % 'guess' size of simulation if it is not provided
            if isempty(opt.N) 
                switch opt.potential_type
                    case "tensor"
                        opt.N = compatible_size(szD(3:end), sza);
                    case "diagonal"
                        opt.N = compatible_size(szD(2:end), sza);
                    case "scalar"
                        opt.N = compatible_size(szD, sza);
                end
            end

            % diffusion simulations are always 4D
            opt.N = extend(opt.N, [1 1 1 1]);

            % verify that dimensions of diffusion coefficient and
            % absorption rate are compatible
            compatible_size(sza, opt.N);
            switch opt.potential_type
                case "tensor"
                    compatible_size(szD(3:end), opt.N);
                    if (szD(1) ~= 3 || szD(2) ~= 3)
                        error("Incorrect dimensions for diffusion tensor field");
                    end
                case "diagonal"
                    compatible_size(szD(2:end), opt.N);
                    if (szD(1) ~= 3)
                        error("Incorrect dimensions for diffusion tensor field");
                    end
                case "scalar"
                    compatible_size(szD, opt.N);
            end

            opt = validate@GridOptions(opt, opt.N);
        end
    end
end
