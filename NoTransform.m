classdef NoTransform
    %NOTRANSFORM Dummy transform class that does not do anything
    %
    properties
        real_signal     % true to indicate that the signal (in r space) is real.
                        % imaginary part is discarded.
    end
    
    methods
        function obj = NoTransform()
        end
        function u = k2r(obj, u, state)  %#ok<INUSD,INUSL>
        end
        function u = r2k(obj, u, state)  %#ok<INUSD,INUSL>
        end
        function u = r2r(obj, u, state)  %#ok<INUSD,INUSL>
        end
    end    
end

