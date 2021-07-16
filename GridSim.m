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
        % required 'options': 
        %   opt.N
        %   opt.unit            e. g. 'mm' or 's'. All dimensions have the same
        %                       unit
        %   opt.pixel_size      size of a pixel (in units). Can be a scalar
        %                       for isotropic pixels, or a vector for 
        %                       anisotropic pixels (i. e. a grid with
        %                       different pitch in different directions)
        %
            defaults.boundaries.periodic = "auto";
            defaults.boundaries.extend = true;
            defaults.boundaries.width = 32;
            defaults.boundaries.filter = @wnd_nutall;
            defaults.potential_type = "scalar"; 
            defaults.crop = true; % otherwise, keeps boundary layers (for debugging)
            opt = set_defaults(defaults, opt);
            obj@AnySim(opt);
            obj.grid = SimGrid(opt.N, opt.boundaries, opt.pixel_size, opt.pixel_unit);
            obj.value_dim = length(N_components);
            value_dims = [N_components, 1, 1];
            obj.N = [value_dims(1:2), obj.grid.N];
            obj.opt.scale_adjuster = @(c, r) obj.adjustScale(c,r);
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
        % Inside the absorbing boundaries, the solution should
        % decay sufficiently to avoid wrap-around artifacts. 
        % Assumining exponential decay, this function returns
        % the minimum decay coefficients needed.
        %
        function mu_min = mu_min(obj)
            mu_min = 10./(obj.grid.boundaries.width .* obj.grid.pixel_size);
            mu_min(obj.grid.boundaries.width == 0) = 0;
            mu_min = max(mu_min, 1E-3 ./ obj.grid.dimensions); % give tiny non-zero minimum value to prevent division by zero in homogeneous media
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
        % Absorbing boundaries are implemented by smoothly tapering
        % G to 0 at the end of the simulation domain. A value of
        % 0 corresponds to maximum absorption. The non-scaled interaction
        % potential at the edges then equals V0+||V||. 
        % However, in the case that the medium is homogeneous and non-
        % absorbing, we have V0=||V||=0, and there would be no
        % absorption at the edges.
        %
        % To fix this problem, we determine the minimum required value 
        % for V0 + ||V|| to ensure sufficient absorption across
        % a boundary with a user-specified width, and include this
        % value as one of the points in V before performing the
        % scaling.
        %
        % mu_min is the minimum absorption coefficient required 
        % in each dimension. Vmax is the current value of
        % V0 + ||V||. If it is already large enough, we just
        % return Vmin = Vmax.
        Vmin = adjustScale(obj, centers, radii)
    end
end
