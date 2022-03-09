classdef TextCallback
    %TEXTCALLBACK Simple callback that just prints a textual progress bar
    % 
    methods
        function obj = TextCallback(opt)
        end
        function obj = prepare(obj, ~)
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

