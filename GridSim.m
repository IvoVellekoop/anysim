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
       mu_min    % Inside the absorbing boundaries, the solution should
                 % decay sufficiently to avoid wrap-around artifacts. 
                 % Assumining exponential decay, this vector returns
                 % the minimum decay coefficients needed in all dimensions.
    end
    methods
        function obj = GridSim(N_components, opt)
        % grid = SimGrid object containing grid spacings and dimensions
        % N_components = length of data vector stored at each grid point
        %                [] for scalar simulations. [N] for vector,
        %                [N,M] for matrix fields (not tested)
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
            defaults.crop = true; % otherwise, keeps boundary layers (for debugging)
            opt = set_defaults(defaults, opt);
            obj@AnySim(opt);
            obj.grid = SimGrid(opt.N, N_components, opt.boundaries, opt.pixel_size, opt.pixel_unit);
            
        
            obj.mu_min = 10./(obj.grid.boundaries.width .* obj.grid.pixel_size);
            obj.mu_min(obj.grid.boundaries.width == 0) = 0;
            obj.mu_min = max(obj.mu_min, 1E-3 ./ obj.grid.dimensions); % give tiny non-zero minimum value to prevent division by zero in homogeneous media
        end
        
        function S = define_source(obj, values, position)
            % SIM.DEFINE_SOURCE(VALUES, POSITION) Defines a source
            % with values specified in the VALUES array, and located
            % at POSITION (in grid points). 
            % For vector simulations, the first element of POSITION
            % specifies at what component of the vector data the source
            % is placed, and the first dimension of VALUES runs along the
            % components
            % 
            % For example, to define an Ez-polarized light source
            % at location 10-12,10,10, in solving 3-D Maxwell's equations,
            % we'd write
            %   sim.define_source(ones(1, 3, 1, 1), [3, 10, 10, 10])
            %
            % Note: the source is converted to an array that is the size of
            %       the full simulation. Theoretically, a better solution
            %       would often be to keep the source as a sparse array,
            %       but MATLAB's has poor support and performace for sparse
            %       arrays at the moment.
            
            if nargin < 3 % defaults to pos = [1,1,1,1...]
                position = [];
            end

            % extend position and size vectors to have length Ndim, and
            % apply offset due to boundaries
            offset = [position - 1, zeros(1, length(obj.grid.N_u) - length(position))]...
                + [zeros(1, length(obj.grid.N_components)) floor(obj.grid.boundaries.width)];
            sz = size(values, 1:length(offset));
            if (any(sz + offset > obj.grid.N_u))
                warning('source does not fit inside simulation window');
            end

            % scale source with matrix Tl, and resize to fit full
            % simulation
            values = padarray(values, offset, 0, 'pre');
            values = padarray(values, obj.grid.N_u - sz - offset, 0, 'post');
            values = fieldmultiply(obj.Tl, values);
            S = data_array(values, obj.opt);
        end
        function c = coordinates(obj, dim)
            % SIM.COORDINATES(D) returns the coordinates in dimension
            % D, cropped to the region of interest. Coordinates are 0 for
            % the start of the ROI, and increasing. The returned array will have
            % size 1 in all dimensions, except in the Dth 
            % dimension, so it is usable for singleton expansion directly
            % (e.g. r = sqrt(sim.coordinates(1).^2 + sim.coordinates(2).^2))
            c = obj.grid.crop(obj.grid.coordinates(dim));
        end

        function u = preconditioner(obj, u)
            % SIM.PRECONDITIONER(U) returns (1-V)(L+1)^(-1)U
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.preconditioner
            u = reshape(u, obj.grid.N_u); 
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
            u = reshape(u, obj.grid.N_u); 
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
            u = reshape(u, obj.grid.N_u); 
            u = operator@AnySim(obj, u);
            u = u(:);
        end
    end
    methods (Access = protected)
        function [u, state] = start(obj)
            u = zero_array(obj.grid.N_u, obj.opt);
            state = State(obj, obj.opt);
        end
        function u = finalize(obj, u, state)  %#ok<INUSD>
            if obj.opt.crop
                u = obj.grid.crop(u, length(obj.grid.N_components));
            end
            u = fieldmultiply(obj.Tr, u);
        end
    end
end
