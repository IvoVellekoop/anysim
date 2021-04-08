classdef Source
    % Object to represent one or more sources for a simulation
    %
    % to create a point source at pos x,y:   s = source([x,y], value);
    % to create a planar source at depth y:  s = source([1,y], ones(N,1));
    % etc.
    %
    % Sources can be added together with the + operator
    % for example:
    %
    % two_points = source([x1,y1], 1) + source([x2, y2], 1);
    %
    %
    % Most functions allow for the use of arrays of source objects
    %
    % Ivo M. Vellekoop 2018
    properties
        %don't access directly, implementation may change!
        ranges   % cell array of index ranges defining volume
        values   % complex amplitude of sources (3-D arrays, dimensions can be 1)
    end
    
    methods
        function obj = Source(values, position, N)
            % Creates a new source at the specified position
            %
            % values = scalar (for a point source), or any array (1-D, 2-D,
            %       3-D, 4-D) describing the amplitude and dimensions
            %       of the source.
            %
            % position = relative position of source within roi of the simulation
            %       Defaults to [1,1,1,1,...] (top-left corner)
            %
            % N = dimension of simulation grid
            %
            
            validateattributes(position, {'numeric'}, {'positive', 'integer'});

            % extend position and size vectors to have length N_dims
            N_dims = length(N);
            position = [position(:); ones(N_dims-length(position), 1)];
            sz = size(values);
            sz = [sz(:); ones(N_dims-length(sz), 1)];
            
            if any(position > N)
                error('The source object is completely outside the simulation grid')
            end
            
            % define index ranges
            last = min(position + sz - 1, N(:));
            ranges = cell(N_dims, 1);
            for d=1:N_dims
                ranges{d} = position(d):last(d); %corresponds to : operator
            end
            obj.ranges = ranges;
            obj.values = values;
        end
                        
        function u_r = apply(obj, u_r, state) %#ok<INUSD>
            % adds source to u_r
            % todo: optimize (for some reason array indexing in Matlab is
            % incredibly slow!!!)
            u_r(obj.ranges{:}) = u_r(obj.ranges{:}) + obj.values;
        end
    end
end

