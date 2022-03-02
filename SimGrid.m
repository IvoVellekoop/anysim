classdef SimGrid
    %SIMGRID Defines coordinate range for a grid-based simulation
    % Ivo M. Vellekoop
    properties
        N       % vector with total number of grid points for all dimensions (_including_ boundaries)
        N_roi   % vector with number of grid points for all dimensions excluding boundaries
        N_dim   % === length(N)
        N_u     % dimensions of data array, includes leading dimensions for components (if N_components not [])
        N_components % size at each grid point
                 % []    = scalar
                 % N     = vector
                 % [N,M] = matrix
        roi_ranges % cell array with index vectors to index the ROI (see crop)
        boundaries % number of grid points used for each boundry 
                   % (can be 1/2 integer number, in this case, on the
                   % left side (low indices) we have floor(boundaries)
                   % and on the right hans side (high indices) we have
                   % ceil(boundaries)
        pixel_size  % vector with step size for all dimensions
        pixel_size_f % vector with step size in Fourier tranformed coordinates
        pixel_unit   % name of the unit of the dimensions (e. g. 'm')
    end
    methods
        function obj = SimGrid(N, N_components, boundaries, pixel_size, pixel_unit)
            % Construct a wave simulation grid object with specified size
            % and boundary options. Automatically adjusts the boundary
            % widths to the closest efficient value.
            %
            % N                   size in pixels _without_ boundaries
            % N_components        size of the data at each grid point ([] =
            %                       scalar, N = vector, [N,M] = matrix)
            % boundaries.width    size in pixels of the added absorbing boundaries 
            %                     (scalar = same for all boudaries,
            %                     vector = possibly different per boundary)
            %                     specify 0 for a periodic boundary.
            % boundaries.extend   when true, boundaries may be made slightly larger
            %                     than specified to allow for an efficient fft.
			%                     Periodic boundaries are never extended.
			% pixel_size    vector containing pixel sizes along the different
            %               dimensions. If there are fewer
            %               entries than dimensions.
            % unit            e. g. 'mm' or 's'. All dimensions have the same
            %                 unit
            
            %% Validate inputs
            validateattributes(N, {'numeric'}, {'positive', 'integer', 'vector'}); 
            validateattributes(boundaries.width, {'numeric'}, {'positive', 'integer'});
            validateattributes(boundaries.extend, {'logical'}, {}); 
            validateattributes(pixel_size, {'numeric', 'vector'}, {'positive'});
            validateattributes(pixel_unit, {'char', 'string'}, {'scalartext'});
            
            % By default, only consider boundaries periodic if the
            % corresponding dimension has size 1.
            if isequal(boundaries.periodic, "auto")
                boundaries.periodic = (N == 1);
            end
            
            % set default values
            if isequal(boundaries.extend, true) % when extend=true, extend all non-periodic boundaries
                boundaries.extend = ~boundaries.periodic;
            end
            if isscalar(boundaries.width) % when scalar boundary width is given, only apply to non-periodic boundaries
                boundaries.width = ~boundaries.periodic * boundaries.width;
            end
            
            %automatically extend to vector of correct size if only a scalar was passed
            obj.N_dim = length(N);
            obj.N_roi = N(:).';
            boundaries.width = SimGrid.extend(boundaries.width(:), obj.N_dim).';
            boundaries.extend = SimGrid.extend(boundaries.extend(:), obj.N_dim).';
            boundaries.periodic = SimGrid.extend(boundaries.periodic(:), obj.N_dim).';
            pixel_size = SimGrid.extend(pixel_size, obj.N_dim).';
            
            % check validity
            if any(boundaries.extend & boundaries.periodic)
                warning('Cannot extend periodic boundaries, "extend" option ignored');
                boundaries.extend = boundaries.extend & ~boundaries.periodic;
            end
            if any(boundaries.periodic & boundaries.width ~= 0)
                warning('Specified boundary width while periodic = true, ignoring periodic');
                boundaries.periodic = boundaries.periodic & boundaries.width == 0;
            end
            if any(~boundaries.periodic & boundaries.width == 0)
                warning('Using a non-periodic boundary without absorbing layer (boundary.width = 0). Results may be unexpected');
            end
            
            %% set up coordinates
            % add boundary size
            obj.boundaries = boundaries;
            obj.N = obj.N_roi + 2 * boundaries.width;
            obj.N(boundaries.extend) = SimGrid.efficient_size(obj.N(boundaries.extend)); %increase size to efficient number for fft
            obj.boundaries.width = (obj.N - obj.N_roi)/2; %boundaries after correction
            
            obj.N_components = N_components;
            obj.N_u = [N_components, obj.N];

            % parse step size and units
            obj.pixel_size = pixel_size;
            obj.pixel_unit = pixel_unit;
            obj.pixel_size_f = 2*pi./(obj.pixel_size .* obj.N);
            
            ranges = cell(obj.N_dim, 1);
            for d=1:obj.N_dim
                ranges{d} = (1:obj.N_roi(d)) + floor(obj.boundaries.width(d));
            end
            obj.roi_ranges = ranges;
        end
        function x = coordinates(obj, d)
            % returns coordinate range for the given dimension
            x = shiftdim(((1:obj.N(d))-floor(obj.boundaries.width(d))-1)*obj.pixel_size(d), 2-d);
        end
        function k = coordinates_f(obj, dimension)
            % returns coordinate range for the given dimension for the
            % Fourier transformed data. Note that the coordinates include a
            % 2pi prefactor (corresponding e.g. to angular frequency omega)
            % and that the coordinates correspond to the data directly after
            % fftn (i.e. without fftshift)
            %
            k = shiftdim(SimGrid.fft_range(obj.N(dimension)) * obj.pixel_size_f(dimension), 2-dimension);
        end
        function X = roi2full(obj, X)
            % Converts a coordinate vector from coordinates with respect to
            % the ROI to coordinates with respect to the full simulation
            % (adding floor(boundaries)).
            % If the coordinate vector is shorter than N_dim, it is
            % right-padded with ones before the conversion
            if length(X) < obj.N_dim
                X(obj.N_dim) = 1;
            end
            X = X + floor(obj.boundaries.width);
        end
        function x = crop(obj, x, d)
            % GRID.CROP(DATA, D) crops the DATA array to the ROI
            % the first D are left untouched (use D=1
            % if the data is a vector field, and 2 if it is a tensor field)
            % This function is compatible with singleton dimensions: 
            % singleton dimensions are kept, not cropped.            
            if nargin < 3
                d = 0;
            end
            validateattributes(d, {'numeric'},{'>=', 0, '<=',2});
            roi = obj.roi_ranges;
            singleton = size(x, (d+1):ndims(x)) == 1;
            roi(singleton) = {1}; % don't crop singleton dimensions
            if d==0
                x = x(roi{:});
            elseif d==1
                x = x(:, roi{:});
            else % d==2
                x = x(:, :, roi{:});
            end
        end
        function M = pad(obj, M, element_dimension, padval)
            % GRID.PAD(DATA, SKIPDIM) appends boundaries to DATA
            % the first SKIPDIM are left untouched (use SKIPDIM=1
            % if the data is a vector field, and 2 if it is a tensor field)
            % PADVAL = padding value or padding method
            % (circular,replicate,symmetric, see padarray)
            if nargin < 3
                element_dimension = 0;
            end
            if nargin < 4
                padval = 'replicate';
            end
            sz = size(M, element_dimension+(1:obj.N_dim));
            if any(sz ~= obj.N_roi & sz ~= 1)
                warning('Padding data that does not have the same size as the simulation window')
            end
            %% Perform padding
            % If padding is needed, first perform singleton expansion for that
            % dimension
            expand = (obj.boundaries.width ~= 0) & (sz == 1);
            M = repmat(M, [ones(1, element_dimension) expand .* (obj.N_roi-1) + 1]); 
            width = [zeros(1, element_dimension) obj.boundaries.width];
            M = padarray(M, floor(width), padval, 'pre');
            M = padarray(M, ceil(width), padval, 'post');

            %% Construct absorption/anti-reflection filter along each of the dimensions
            for d=1:obj.N_dim
                w = obj.boundaries.width(d);
                if w>0 % construct vector with boundaries on both sides, and 1 in between
                    left_boundary = obj.boundaries.filter(floor(w));
                    right_boundary = obj.boundaries.filter(ceil(w));
                    full_filter = [left_boundary(:); ones(obj.N(d)-floor(w)-ceil(w), 1); flipud(right_boundary(:))];

                    % transpose vector to proper dimension
                    full_filter = shiftdim(full_filter, 1-d-element_dimension);
                    
                    % apply filter along this dimension
                    M = M .* full_filter;
                end
            end
        end
            
        function d = dimensions(obj)
            % GRID.DIMENSIONS returns a vector with the total size of the
            % simulation grid (in units)
            d = obj.N .* obj.pixel_size;
        end
    end
	methods(Static)
        function range = symrange(N)
            % Cnstructs a range from -floor(N/2) to ceil(N/2)-1
            % This is the natural way to define coordinates in
            % MATLAB:
            % - element 0 is always part of the range
            % - the step size is exacty 1/N
            % - the range is compatible with fft. 
            %
            % When N is odd, the range is symmetric around 0 (from -floor(N/2) to floor(N/2)).
            % When N is even, the range is almost symmetric and includes 0
            % (from -N/2 to N/2-1)
            %
            range = -floor(N/2):ceil(N/2)-1;
        end
        function range = fft_range(N)
            % constructs a range of numbers that is compatible
            % with the coordinates after a fft.
            % implemented as ifftshift(symrange(N))
            range = ifftshift(SimGrid.symrange(N));
        end
        function sz = efficient_size(min_size)
            % returns nearest size greater than or equal to min_size
            % for which the fft is efficient.
            % cuFFT is efficient for sizes that can be factored as 2^a *
            % 3^b * 5^c * 7^d * 11^e (tested empirically using test_efficient_size)
            %
            sz = min_size;
            for s_i=1:length(sz)
                s = sz(s_i);
                f = factor(s);
                while f(end) > 11 || (length(f)>2 && f(end-1) > 5)
                    s = s + 1;
                    f = factor(s);
                end
                sz(s_i) = s;
            end
        end
        function test_efficient_size()
            s_range = 1:200;
            times = zeros(numel(s_range), 1);
            efficient = zeros(numel(s_range), 1, 'logical');
            for s_i=1:numel(s_range)
                s = s_range(s_i);
                f = factor(s);
                efficient(s_i) = true;
                if f(end) > 11 || (length(f)>2 && f(end-1) > 5)
                    efficient(s_i) = false;
                end
                data = gpuArray((1+1i)*ones([s,s,s], 'single'));
                times(s_i) = gputimeit(@() fftn(data)) / (s * log2(s));
                plot(s_range(efficient), times(efficient), '.'); 
                hold on;
                plot(s_range(~efficient), times(~efficient), '*'); 
                hold off;
                drawnow;
            end
        end
        function L = fix_edges_hermitian(L, dimensions)
            % For dimensions with even sizes, there is one coordinate
            % for which there is no mirror point (i.e. -x exists
            % but x is not in the coordinate range)
            % When trying to make a Hermitian signal this single
            % data point breaks the Hermitian symmetry.
            % To solve this problem, we replace the point at
            % -x by its real part: L[-x] = Re(L[x]) so that the
            % resulting signal has Hermitian symmetry.
            % This function operates along the specified 'dimensions'
            % and only has an effect when the size(L, dimension) is even
            for d=dimensions
                N = size(L, d);
                if mod(N, 2) == 0
                    % construct indexing operation L(:,:,end/2+1,:)
                    ind.type = '()';
                    [ind.subs{1:ndims(L)}] = deal(':');
                    ind.subs{d} = N/2 + 1;
                    L = subsasgn(L, ind, real(subsref(L, ind)));
                %    L = subsasgn(L, ind, 0);
                end
            end
        end
    end
    methods(Static, Access=private)
        function x = extend(x, L)
            % repeats the last element in the vector until length L is
            % reached
            Lorig = size(x, 1);
            for l = L:-1:(Lorig+1)
                x(l,:) = x(Lorig,:);
            end
        end
    end
end



