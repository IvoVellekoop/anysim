classdef PrintIterationCallback
    %PRINTITERATIONCALLBACK Simple callback that just prints the iteration number
    % 
    methods
        function obj = PrintIterationCallback(sim, opt)
        end
        
        function call(obj, u, state)
            disp(state.iteration)
        end
    end
end

