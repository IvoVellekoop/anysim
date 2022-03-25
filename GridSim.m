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
        function obj = GridSim(N, opt)
            arguments
                % Dimensions of the grid in voxels. When empty, this should be 
                % determined automatically from the passed potential array.
                N (1,:) {mustBePositive}
                opt (1,1) {mustBeA(opt, "AnySimOptions"), mustBeA(opt, "GridOptions")}
            end

            obj@AnySim(opt);
            obj.grid = Grid(N, opt);
            
            obj.mu_min = 10./(obj.grid.boundaries_width .* obj.grid.pixel_size);
            obj.mu_min(obj.grid.boundaries_width == 0) = 0;
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
            valuedim = length(obj.grid.N_components);
            offset = [position - 1, zeros(1, length(obj.grid.N_u) - length(position))]...
                + [zeros(1, valuedim) floor(obj.grid.boundaries_width)];
            sz = size(values, 1:length(offset));
            does_not_fit = sz + offset > obj.grid.N_u;
            if any(does_not_fit(1:valuedim))
                error('source does not fit inside simulation window, are you missing a leading singleton dimension? (shiftdim -1)');
            end
            if any(does_not_fit)
                warning('source does not fit inside simulation window');
            end

            % scale source with matrix Tl, and resize to fit full
            % simulation
            values = padarray(values, offset, 0, 'pre');
            values = padarray(values, obj.grid.N_u - sz - offset, 0, 'post');
            values = fieldmultiply(obj.Tl, values);
            S = obj.data_array(values);
        end
        function D = data_array(obj, data)
            %DATA_ARRAY Converts an array to the proper data type (double/single, gpu or not)
            %   DATA_ARRAY(X, OPT)      converts data in X to formatting specified in OPT
            %      OPT.precision   = 'single' or 'double'
            %      OPT.gpu_enabled = 'true' or 'false'
            %
            switch obj.opt.precision
            case 'single'
                D = single(data);
            case 'double'
                D = double(data);
            otherwise
                error('Precision should be single or double. %s is not supported', opt.precision);
            end
            if obj.opt.gpu_enabled
                D = gpuArray(D);
            else
                D = gather(D);
            end
        end
        function D = zero_array(obj, N)
            %ZERO_ARRAY Constructs a zero array of size N with proper formatting (double/single, gpu or not)
            %   OPT.precision   = 'single' or 'double'
            %   OPT.gpu_enabled = 'true' or 'false'
            %
            Dscalar = obj.data_array(0); % avoid code duplication with data_array
            D = zeros([N 1], 'like', Dscalar);
        end
        function u = preconditioner(obj, u)
            % SIM.PRECONDITIONER(U) returns (1-V)(L+1)^(-1)U
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.preconditioner
            u = reshape(u, [obj.grid.N_u, 1]); 
            u = preconditioner@AnySim(obj, u);
            u = u(:);
        end
        function [f, state] = preconditioned(obj)
            % SIM.PRECONDITIONED
            %
            % See AnySim.preconditioned
            % Returns a function that computes the preconditioned operator
            % A
            % Also see AnySim.preconditioned
            [A, state] = preconditioned@AnySim(obj);
            f = @(u, varargin) reshape(A(reshape(u, [obj.grid.N_u, 1]), varargin{:}), [], 1); 
        end
        function [f, state] = operator(obj)
            % SIM.OPERATOR(U) Returns (L+V)U
            %
            % For compatibility with MATLAB built in algorithms
            % such as GMRES, the input and output are column vectors
            %
            % Also see AnySim.operator
            [A, state] = operator@AnySim(obj);
            f = @(u, varargin) reshape(A(reshape(u, [obj.grid.N_u, 1]), varargin{:}), [], 1); 
        end
        function u = finalize(obj, u)
            u = reshape(u, obj.grid.N_u);
            u = obj.grid.crop(u);
            u = fieldmultiply(obj.Tr, u);
        end
    end
    methods (Access = protected)
        function [u, state] = start(obj)
            u = obj.zero_array(obj.grid.N_u);
            state = State(obj, obj.opt);
        end
    end
end
