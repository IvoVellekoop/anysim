classdef HelmholtzSimOptions < AnySimOptions
    %HELMHOLTZSIMOPTIONS Options for simulation of the Helmholtz equation
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % Dimensions of the grid in voxels, excluding boundaries.
        % When empty, the size is determined automatically from n.
        N (1,:) {mustBePositive} = []
        wavelength (1,1) double = 1
        grid GridOptions = GridOptions(pixel_size = 0.25, pixel_unit = 'λ')
    end
    methods 
        function obj = HelmholtzSimOptions(opt)
            arguments 
                opt.?HelmholtzSimOptions
            end
            obj.alpha = 0.75; % empirical default value for wave operators
            obj = copy_properties(obj, opt);
        end
        function opt = validate(opt, szn)
            % the simulation is always scalar
            opt.grid.N_components = [];

            % 'guess' size of simulation if it is not provided
            if isempty(opt.N) 
                opt.N = szn;
            end
        end
    end
end
