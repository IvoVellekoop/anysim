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
       grid      % simulation grid
       N         % dimensions of data array + 2 leading dimensions for components (if N_components > 0)
       value_dim % dimensionality of values stored at each grid point
                 % 0 = scalar
                 % 1 = vector
                 % 2 = matrix                 
    end
    methods
        function obj = GridSim(N_components, opt)
        % grid = SimGrid object containing grid spacings and dimensions
        % N_components = length of data vector stored at each grid point
        %                [] for scalar simulations. [N] for vector,
        %                [N,M] for matrix.
        % required options: 
        %   opt.N
        %   opt.pixel_size
        %
            defaults.boundaries.periodic = "auto";
            defaults.boundaries.extend = true;
            defaults.boundaries.width = 32;
            defaults.boundaries.filter = @wnd_nutall;
            defaults.potential_type = "scalar"; 
            defaults.crop = true; % otherwise, keeps boundary layers (for debugging)
            opt = set_defaults(defaults, opt);
            obj@AnySim(opt);
            obj.grid = SimGrid(opt.N, opt.boundaries, opt.pixel_size);
            obj.value_dim = length(N_components);
            value_dims = [N_components, 1, 1];
            obj.N = [value_dims(1:2), obj.grid.N];
        end
        
        function S = define_source(obj, values, position)
            % SIM.DEFINE_SOURCE(VALUES, POSITION) Defines a source
            % with values specified in the VALUES array, and located
            % at POSITION (in grid points). 
            % For vector simulations, the first element of POSITION
            % specifies at what component of the vector data the source
            % is placed. For inherently scalar simulations (e.g.
            % Helmholtz), the POSITION vector just specifies coordinates
            % (e.g. x,y,z,t).
            % 
            % For example, to define an Ez-polarized light source
            % at location 10,10,10, in solving 3-D Maxwell's equations,
            % we'd write
            %   sim.define_source(1, [3, 10, 10, 10])
            %
            % Effectively, the data is stored at location
            % source(position(1) + (1:size(values,1)),
            %        position(2) + (1:size(values,2)),
            %        position(3) + (1:size(values,3)), etc.) = values
            %
            if nargin < 3 % defaults to pos = [1,1,1,1...]
                position = ones(1, length(obj.N));
            else
                % convert to internal representation coordinates
                % for scalar simulations: prepend 1,1
                %   change from [x, y, ...] to [1, 1, x, y , ...]
                %
                % for vector simulations: 
                %   change from [v, x, y, ...] to [v, 1, x, y , ...]
                %
                if obj.value_dim == 1
                    position = [position(1), 1, position(2:end)];
                else % value_dim == 0
                    position = [1, 1, position];
                end
            end
            values = obj.to_internal(values, obj.value_dim);
            
            % scale source with matrix Tl. If only values for
            % part of the components were specified, make sure
            % to use the correct part of the matrix!
            if ~isdiag(obj.medium.Tl) || obj.value_dim == 2
                error('only works for diagonal matrices at the moment!');
            end
            source_components = position(1)-1+(1:size(values,1));
            if source_components(end) > obj.N(1)
                error("Specified more source components than possible");
            end
            scale = obj.medium.Tl(source_components, source_components);
            values = pagemtimes(scale, values);
                
            % translates the source to account for the boundaries in the grid
            % gives error if the source is completely outside the grid
            % todo: warn if source overlaps boundaries?
            position = [position(1:2) obj.grid.roi2full(position(3:end))];
            S = Source(values, position, obj.N);
        end
        function u = to_internal(obj, u, dim)
            % convert u from external to internal representation
            % scalar:   internal:   [1  1  Nx Ny ...]
            %           external:   [Nx Ny ...]
            % vector:   internal:   [Nc 1  Nx Ny ...]
            %           external:   [Nc Nx Ny ...]
            % matrix:   internal:   [Nn Nm Nx Ny ...]
            %           external:   [Nn Nm Nx Ny ...]
            if nargin < 3
                dim = obj.value_dim;
            end
            if dim == 0
                u = shiftdim(u, -2);
            elseif dim == 1
                sz = size(u);
                u = reshape(u, [sz(1), 1, sz(2:end)]);
            %elseif dim == 2: no format change needed
            end                
        end
        function u = to_external(obj, u, dim)
            % convert u from external to internal representation
            % scalar:   internal:   [1  1  Nx Ny ...]
            %           external:   [Nx Ny ...]
            % vector:   internal:   [Nc 1  Nx Ny ...]
            %           external:   [Nc Nx Ny ...]
            % matrix:   internal:   [Nn Nm Nx Ny ...]
            %           external:   [Nn Nm Nx Ny ...]
            if nargin < 3
                dim = obj.value_dim;
            end
            if dim == 0
                u = shiftdim(u, 2);
            elseif dim == 1
                sz = size(u);
                u = reshape(u, [sz(1), sz(3:end)]);
            %elsif dim == 2: no format change needed
            end                
        end
        function u = preconditioner(obj, u)
            % SIM.PRECONDITIONER(U) returns (1-V)(L+1)^(-1)U
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.preconditioner
            u = reshape(u, obj.N); 
            u = preconditioner@AnySim(obj, u);
            u = u(:);
        end
        function u = preconditioned(obj, u)
            % SIM.PRECONDITIONED(U)
            %
            % See AnySim.preconditioned
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.preconditioner
            u = reshape(u, obj.N); 
            u = preconditioned@AnySim(obj, u);
            u = u(:);
        end
        function u = operator(obj, u)
            % SIM.OPERATOR(U) Returns (L+V)U
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.operator
            u = reshape(u, obj.N); 
            u = operator@AnySim(obj, u);
            u = u(:);
        end
    end
    
    
    methods (Access = protected)
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
            % convert to internal representation (with trailing 1
            % dimensions for scalar)
            if obj.opt.potential_type == "scalar" %convert to diagonal matrix
                dim = 0;
            elseif obj.opt.potential_type == "diagonal" %convert to diagonal matrix
                dim = 1;
            elseif obj.opt.potential_type == "tensor"
                dim = 2;
            end                
            
            Vraw = obj.to_internal(Vraw, dim);
            Vmax = max(Vraw, [], 3:max(ndims(Vraw), 3));
            Vmin = obj.analyzeDimensions(Vmax);
                
            % compute the current maximum value of V (before scaling),
            if obj.opt.potential_type == "tensor"
                validateattributes(Vraw, {'numeric'}, {'nrows', Nc, 'ncols', Nc});
                medium = TensorMedium(Vraw, Vmin, obj.grid, obj.opt);
            elseif obj.opt.potential_type == "diagonal"
                medium = DiagonalMedium(Vraw, Vmin, obj.grid, obj.opt);
            elseif obj.opt.potential_type == "scalar" %convert to diagonal matrix
                medium = ScalarMedium(Vraw, max(Vmin), obj.grid, obj.opt);
            end
        end
        function [u, state] = start(obj)
            u = zero_array(obj.N, obj.opt);
            state = State(obj, obj.opt);
        end
        function u = finalize(obj, u, state)  %#ok<INUSD>
            if obj.opt.crop
                u = obj.grid.crop(u, 2);
            end
            u = pagemtimes(obj.medium.Tr, u);
            u = obj.to_external(u, obj.value_dim); %remove spurious dimensions
        end
    end
    methods (Abstract, Access = protected)
        Vmin = analyzeDimensions(obj, Vmax)
        % Analyzes the dimensions and spacings of the simulation
        % grid. Issues a warning if the grid spacing is larger
        % than the size of the smallest features.
        % Also computes a minimum value for the scattering potential
        % to ensure that the Green's function does not extend beyond
        % the size of the simulation window.
    end
end
