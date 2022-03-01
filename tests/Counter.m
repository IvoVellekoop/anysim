classdef Counter < handle
    %% Helper class to count the number of times an anynymous function is invoked
    % Usage:
    % c = Counter
    % f = @(input) c.inc(some_function(input))
    % ... (use f)...
    % disp(c.i);  % show number of times f was called
    % c.reset
    %
    % inc is a function that passes the input argument
    % unchanged and increments the counter.
    properties 
        i
    end
    methods
        function obj = Counter()
            obj.i = 0;
        end
        function dummy = inc(obj, dummy)
            obj.i = obj.i + 1;
        end
        function reset(obj)
            obj.i = 0;
        end
    end
end