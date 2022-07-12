classdef PantographOptions < AnySimOptions & GridOptions
    %PANTOGRAPHOPTIONS Options for solving the pantograph equation
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % Dimensions of the grid in voxels, excluding boundaries.
        % When empty, the size is determined automatically from D and a.
        N (1,:) {mustBePositive} = []

        % Set to false to solve a non-accretive system. This solves
        % the normal equations (anti-symmetrized system) instead.
        accretive (1, 1) logical = true
    end
    methods 
        function obj = PantographOptions(opt)
            arguments 
                opt.?PantographOptions
            end
            %% Set defaults
            obj.pixel_unit = 's';
            opt.gpu_enabled = false;
            % disable gpu. For 'Pantograph.convolve' is is very inefficient, and
            % for PantographF we run into MATLAB bugs (case 05553718)
            % involving sparse gpu matrices.
            %
            opt.precision = 'double'; % needed because we work with sparse matrices in PantographF
            obj = copy_properties(obj, opt);
        end
        function opt = validate(opt, sza, szb)
            % the pantograph is always a 1-D equation with scalar components
            if sza(2) > 1 || szb(2) > 1 || length(sza) > 2 || length(szb) > 2 || length(opt.N)>1
                error('alpha and beta must be 1-D arrays or scalars');
            end
            opt.N_components = [];

            % 'guess' size of simulation if it is not provided
            if isempty(opt.N) 
                opt.N = compatible_size(sza, szb);
            end

            % verify that dimensions of alpha and beta are compatible
            compatible_size(sza, opt.N);
            compatible_size(szb, opt.N);
            opt = validate@GridOptions(opt, opt.N);
        end
    end
end
