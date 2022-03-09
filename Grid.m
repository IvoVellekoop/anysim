classdef Grid < GridOptions
    %SIMGRID Defines coordinate range for a grid-based simulation
    % Ivo M. Vellekoop
    properties
        N       % size of the simulation grid in voxels including absorbing boundaries and padding
        N_roi   % vector with number of grid points for all dimensions excluding boundaries
        N_dim   % === length(N)
        N_u     % dimensions of data array, includes leading dimensions for components (if N_components not [])
        value_dim % === length(N_components)
        roi_ranges % cell array with index vectors to index the ROI (see crop)
        pixel_size_f % pixel size after Fourier transform
    end
    methods
        function obj = Grid(N, opt)
            arguments
                % size of the simulation grid in voxels,
                % excluding absorbing boundaries or padding
                N (1,:) {mustBeInteger, mustBePositive}
                opt (1,1) GridOptions
            end

            % copy options into this object
            obj = copy_properties(obj, opt);
            obj.N_dim = length(N);
            obj.N_roi = N;
            obj.value_dim = length(opt.N_components);

            %% set up coordinates
            % add boundary size
            % increase size to efficient number for fft
            obj.N = obj.N_roi + round(2 * obj.boundaries_width);
            obj.N(obj.boundaries_extend) = Grid.efficient_size(obj.N(obj.boundaries_extend));
            obj.boundaries_width = (obj.N - obj.N_roi)/2; % boundaries after correction
            obj.N_u = [obj.N_components, obj.N];
            obj.pixel_size_f = 2*pi./(obj.pixel_size .* obj.N);

            % compute ranges for selecting ROI. Todo: create subsref
            % structure instead?
            ranges = cell(obj.N_dim, 1);
            for d=1:obj.N_dim
                ranges{d} = (1:obj.N_roi(d)) + floor(obj.boundaries_width(d));
            end
            obj.roi_ranges = ranges;
        end
        function x = coordinates(obj, d, edges)
            arguments
                obj
                d (1,1) { mustBeInteger, mustBePositive }
                edges { mustBeMember(edges, ["full", "crop", "auto", "nan"])} = "auto"
            end
            % returns coordinate range for the given dimension. If
            % the crop_to_roi option is set to true, the coordinates
            % are cropped to the roi, unless full = true
            % full: include boundaries
            % crop: remove boundaries
            % auto: crop_to_roi ? "crop" : "full"
            % nan: put 'nan' in boundaries
            %
            x = shiftdim(((1:obj.N(d))-floor(obj.boundaries_width(d))-1)*obj.pixel_size(d), 2-d);
            switch edges
                case "crop"
                    x = x(obj.roi_ranges{d});
                case "auto"
                    if obj.crop_to_roi
                        x = x(obj.roi_ranges{d});
                    end
                case "nan"
                    nans = zeros(size(x)) + nan;
                    nans(obj.roi_ranges{d}) = 0;
                    x = x + nans;
            end
        end
        
        function k = coordinates_f(obj, dimension)
            % returns coordinate range for the given dimension for the
            % Fourier transformed data. Note that the coordinates include a
            % 2pi prefactor (corresponding e.g. to angular frequency omega)
            % and that the coordinates correspond to the data directly after
            % fftn (i.e. without fftshift)
            %
            k = shiftdim(Grid.fft_range(obj.N(dimension)) * obj.pixel_size_f(dimension), 2-dimension);
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
            X = X + floor(obj.boundaries_width);
        end
        function x = crop(obj, x)
            % GRID.CROP(DATA, D) crops the DATA array to the ROI
            % Note: the first obj.N_dim dimensions are left untouched,
            % Note: singleton dimensions are kept, not cropped.
            % Note: if crop_to_roi = false, this function just returns x
            if ~obj.crop_to_roi
                return;
            end
            d = obj.value_dim;
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
            expand = (obj.boundaries_width ~= 0) & (sz == 1);
            M = repmat(M, [ones(1, element_dimension) expand .* (obj.N_roi-1) + 1, 1]); 
            width = [zeros(1, element_dimension) obj.boundaries_width];
            M = padarray(M, floor(width), padval, 'pre');
            M = padarray(M, ceil(width), padval, 'post');

            %% Construct absorption/anti-reflection filter along each of the dimensions
            for d=1:obj.N_dim
                w = obj.boundaries_width(d);
                if w>0 % construct vector with boundaries on both sides, and 1 in between
                    left_boundary = obj.boundaries_window(floor(w));
                    right_boundary = obj.boundaries_window(ceil(w));
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
            range = ifftshift(Grid.symrange(N));
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



