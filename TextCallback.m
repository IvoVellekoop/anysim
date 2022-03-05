classdef TextCallback
    %TEXTCALLBACK Simple callback that just prints a textual progress bar
    % 
    methods
        function obj = TextCallback(sim, opt)
        end
        
        function call(obj, u, r, state)
            if state.iteration == 1
                fprintf('\n');
            else
                fprintf('*');
            end
        end
    end
end

