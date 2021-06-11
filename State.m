classdef State < dynamicprops
    %STATE Holds state variables (variables that change while running the
    %   algorithm.)
    %
    %   (c) 2019. Ivo Vellekoop    
    properties
        iteration;
        start_time;
        termination_condition_interval;
        termination_condition;
        callback_interval;
        callback;
        running;
        diffs;
        diff_its; % iteration numbers at which diffs were reported
    end
    
    methods
        function obj = State(sim, opt)
            obj.iteration = 1; %current iteration number
            obj.start_time = cputime; %only measure time actually used by MATLAB
            obj.running = true;
            obj.callback_interval = opt.callback.interval;
            if ~isempty(opt.callback.handle)
                obj.callback = opt.callback.handle(sim, opt.callback);
            end
            obj.diffs = [];
            obj.termination_condition_interval = opt.termination_condition.interval; % how often to analyze the update du (for termination condition and/or visual feedback)
            obj.termination_condition = opt.termination_condition.handle(sim, opt.termination_condition);
        end
        function next(obj, u)
            if mod(obj.iteration-1, obj.callback_interval) == 0 && ~isempty(obj.callback)
            %    fprintf("Iteration: %d, total time: %g\n", obj.iteration, cputime-obj.start_time);
                obj.callback.call(u, obj);
            end
            obj.iteration = obj.iteration + 1;
        end
        function needed = needs_report(obj)
            % STATE.NEEDS_REPORT returns true if during this iteration
            % additional information should be reported (currently only
            % the magnitude of du, see report_diff)
            %
            needed = mod(obj.iteration-1, obj.termination_condition_interval) == 0;
        end
        function report_diff(obj, diff)
            % STATE.REPORT_DIFF(DIFF) 
            % This callback is called from the medium.apply function
            % when du = u^(k) - u^(k-1) is computed.
            obj.diffs = [obj.diffs, diff];
            obj.diff_its = [obj.diff_its, obj.iteration];
            obj.running = ~obj.termination_condition.call(obj);
        end
    end
    
end

