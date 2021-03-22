classdef GridSim < AnySim
    %GRIDSIM Base class for all grid-based simulations
    %   Base class for all simulations that are performed on a 
    %   regular grid.
    %   For simplicity, grids are always 4-dimensional: (x,y,z,t)
    %   Since each grid point may hold a vector of values, the
    %   data array is 5-dimensional.
    %   If a dimension is not used, it has size 1 
    %   For example, a stationary scalar 3-D simulation has
    %   a data array of size [1, Lx, Ly, Lz, 1]
    %
    %   The grid spacings are stored in the 'grid' object
    %
    %   (c) 2019. Ivo Vellekoop
    properties
       grid     % simulation grid
       N        % dimensions of data array
    end
    methods
        function obj = GridSim(N_components, opt)
        % grid = SimGrid object containing grid spacings and dimensions
        % N_components = length of data vector stored at each grid point
        % opt = options
            opt = set_defaults(GridSim.defaults, opt);
            obj@AnySim(opt);
            obj.grid = SimGrid(opt.N, opt.boundaries, opt.pixel_size);
            obj.N = [N_components obj.grid.N];
        end
        
        % Constructs the medium operator from the potential matrix Vraw
        % The meaning of Vraw depends on opt.potential_type:
        %   "scalar", Vraw has size N (after singleton dimension expansion)
        %   "diagonal", Vraw has size N_components x N. The components are
        %               elements of a diagonal matrix.
        %   "tensor", Vraw has size N_components x N_components x N.
        %               So, for each grid point we have a full matrix.
        %
        function medium = makeMedium(obj, Vraw)
            Nc = obj.N(1);
            % compute the current maximum value of V (before scaling),
            if obj.opt.potential_type == "tensor"
                Vmax = max(Vraw, [], 3:ndims(Vraw));
                Vmin = obj.analyzeDimensions(Vmax);
                validateattributes(Vraw, {'numeric'}, {'nrows', Nc, 'ncols', Nc});
                medium = TensorMedium(Vraw, Vmin, obj.grid, obj.opt);
            elseif obj.opt.potential_type == "diagonal"
                Vmax = max(Vraw, [], 2:ndims(Vraw));
                Vmin = obj.analyzeDimensions(Vmax);
                validateattributes(Vraw, {'numeric'}, {'nrows', Nc});
                medium = DiagonalMedium(Vraw, Vmin, obj.grid, obj.opt);
            elseif obj.opt.potential_type == "scalar" %convert to diagonal matrix
                Vmax = max(Vraw(:));
                Vmin = obj.analyzeDimensions(Vmax);
                validateattributes(Vraw, {'numeric'});
                medium = ScalarMedium(Vraw, max(Vmin), obj.grid, obj.opt);
            end
        end
        
        function S = define_source(obj, values, position)
            % SIM.DEFINE_SOURCE(VALUES, POSITION) Defines a source
            % with values specified in the VALUES array, and located
            % at POSITION (in grid points). 
            % The first element of POSITION specifies at what
            % component of the vector data the source is placed.
            % For example, to define an Ez-polarized light source
            % at location 10,10,10, in solving 3-D Maxwell's equations,
            % we'd write
            %   sim.define_source(1, [3, 10, 10, 10])
            % For scalar simulations, the first element should be 1.
            %
            % Effectively, the data is stored at location
            % source(position(1) + (1:size(values,1)),
            %        position(2) + (1:size(values,2)),
            %        position(3) + (1:size(values,3)), etc.) = values
            %
            
            % scales the source values by -Tl
            validateattributes(position, {'numeric'}, {'positive', 'integer'});
            source_components = position(1)-1 + (1:size(values, 1)); 
            if source_components(end) > obj.N(1)
                error('Invalid position for source. Note: dimension should be N_components * Nx * Ny * Nz * Nt');
            end
            
            % compute Tl * source
            if ~isdiag(obj.medium.Tl)
                error('only works for diagonal matrices!');
            end
            scale = obj.medium.Tl(source_components, source_components);
            values = pagemtimes(scale, values);
            
            % translates the source to account for the boundaries in the grid
            % gives error if the source is completely outside the grid
            % todo: warn if source overlaps boundaries?
            pos = [position(1) obj.grid.roi2full(position(2:end))];
            
            S = Source(values, pos, obj.N);
        end
    end
    methods (Access=protected)
        function [u, state] = start(obj)
            u = zero_array(obj.N, obj.opt);
            state = State(obj, obj.opt);
        end
        function u = finalize(obj, u, state)  %#ok<INUSD>
            u = obj.grid.crop(u, 1);
            u = pagemtimes(obj.medium.Tr, u);
        end
    end
    methods (Static)
        function opt = defaults
            opt = AnySim.defaults;
            opt.boundaries.periodic = "auto";
            opt.boundaries.extend = true;
            opt.boundaries.width = 32;
            opt.boundaries.filter = @wnd_nutall;
            opt.potential_type = "scalar";         
        end
    end
end
