classdef HelmholtzSimOptions < AnySimOptions & GridOptions
    %HELMHOLTZSIMOPTIONS Options for simulation of the Helmholtz equation
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % Dimensions of the grid in voxels, excluding boundaries.
        % When empty, the size is determined automatically from n.
        N (1,:) {mustBePositive} = []
        wavelength (1,1) double = 1
        legacy_mode (1,1) logical = false % for comparison with original wavesim algorithm
    end
    methods 
        function obj = HelmholtzSimOptions(opt)
            arguments 
                opt.?HelmholtzSimOptions
            end
            % override some default values:
            obj.alpha = 0.75; % empirical default value for wave operators
            obj.pixel_size = 0.25;
            obj.pixel_unit = 'λ';
            obj = copy_properties(obj, opt);
        end
        function opt = validate(opt, szn)
            % the simulation is always scalar
            opt.N_components = [];

            % 'guess' size of simulation if it is not provided
            if isempty(opt.N) 
                opt.N = szn;
            else
                % if a size is specified, check if the refractive index
                % array has a compatible size
                sz_len = max(length(szn), length(opt.N));
                szpad = [szn(:); ones(1, sz_len - length(szn))];
                Npad = [opt.N(:); ones(1, sz_len - length(opt.N))];
                if any(Npad ~= szpad & szpad ~= 1)
                    error("opt.N and size of refractive index array don't match");
                end
            end
            opt = validate@GridOptions(opt, opt.N);
        end
    end
end
