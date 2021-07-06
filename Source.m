classdef Source
    % Object to represent one or more sources for a simulation
    %
    % to create a point source at pos x,y:   s = source([x,y], value);
    % to create a planar source at depth y:  s = source([1,y], ones(N,1));
    % etc.
    %
    % todo: Sources can be added together with the + operator
    % for example:
    % two_points = source([x1,y1], 1) + source([x2, y2], 1);
    % todo: allow for the use of arrays of source objects
    %
    % Ivo M. Vellekoop 2021
    properties
        %don't access directly, implementation may change!
        ranges   % cell array of index ranges defining volume
        values   % complex amplitude of sources (3-D arrays, dimensions can be 1)
        N        % dimension of the simulation grid
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
            %       the position vector should have the same length as N
            % N = dimension of simulation grid in internal representation
            %
            
            validateattributes(position, {'numeric'}, {'positive', 'integer'});

            obj.N = N(:).';
            if any(position > obj.N)
                error('The source object is completely outside the simulation grid')
            end
            
            % define index ranges
            N_dims = length(obj.N);
            sz = size(values, 1:N_dims); % extend size vector to have length N_dims
            last = min(position + sz - 1, obj.N);
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
        
        function b = to_array(obj)
            % returns the source object as a full array
            b = zeros(obj.N, 'like', obj.values);
            b = obj.apply(b, []);
        end
    end
end

