classdef GridOptions
    %GRIDOPTIONS Options for describing a simulation grid and its
    %boundaries
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % Vector with units foreach dimensions "mm" or "s".
        % Missing values: repeats last value for remaining dimensions
        pixel_unit (1,:) string {mustBeNonempty} = "-"

        % Size of a voxel (in units).
        % Missing values: repeats last value for remaining dimensions
        pixel_size (1,:) {mustBePositive, mustBeNonempty} = 1

        % size of data vector stored at each grid point
        % [] for scalar simulations. [N] for vector,
        % [N,M] for matrix fields (not tested)
        N_components {mustBeInteger, mustBePositive} = []
        
        % total size in pixels of the added absorbing boundaries for each dimension 
        % the boundary width is added to both sides.
        % Missing values: repeat last value (take as 32 when missing),
        % but set to 0 for singleton dimensions (these are treated as periodic).
        boundaries_width (1,:) {mustBeNonnegative} = []
        
        % when true, boundaries may be made slightly larger than specified
        % to allow for an efficient fft. Periodic boundaries are never extended.
        % Missing values: defaults to 'true' for dimensions with a non-zero boundary_width.
		boundaries_extend logical = []

        % function that is used to compute the absorbing boundary layer
        boundaries_window = @wnd_nutall;

        % set to true (default) to remove boundary and padding layers
        % for instance in a call to coordinates
        crop_to_roi (1,1) logical = true
    end

    methods 
        function obj = GridOptions(opt)
            arguments 
                opt.?GridOptions
            end
            obj = copy_properties(obj, opt);
        end
        
        function opt = validate(opt, N)
            % validates options and sets defaults for given size vector N
    
            % A boundary is considered periodic if a boundary_width of 0 is
            % specified. When no boundary width is specified, a dimension
            % is considered periodic if it's extent is 1 voxel.
            opt.pixel_unit = extend(opt.pixel_unit, repmat(opt.pixel_unit(end), 1, length(N)));
            opt.pixel_size = extend(opt.pixel_size, repmat(opt.pixel_size(end), 1, length(N)));
            if ~isempty(opt.boundaries_width)
                default_width = opt.boundaries_width(end);
            else
                default_width = 32;
            end
            opt.boundaries_width = extend(opt.boundaries_width, (N~=1) * default_width);
            opt.boundaries_extend = extend(opt.boundaries_extend, opt.boundaries_width ~= 0);
            
            % check validity
            if any(opt.boundaries_width == 0 & opt.boundaries_extend)
                warning('Applying zero-padding in a direction where the simulation is periodic, is this desired?')
            end            
        end
    end    
end
