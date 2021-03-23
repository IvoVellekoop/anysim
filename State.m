classdef State < dynamicprops
    %STATE Holds state variables (variables that change while running the
    %   algorithm.)
    %
    %   (c) 2019. Ivo Vellekoop    
    properties
        iteration;
        start_time;
        termination_condition;
        callback_interval;
        callback;
        running;
        diffs;
    end
    
    methods
        function obj = State(sim, opt)
            obj.iteration = 1; %current iteration number
            obj.start_time = cputime; %only measure time actually used by MATLAB
            obj.termination_condition = opt.termination_condition;
            obj.running = true;
            obj.callback_interval = opt.callback.interval;
            obj.callback = opt.callback.handle(sim, opt.callback);
            obj.diffs = [];
        end
        function next(obj, u)
            if mod(obj.iteration, obj.callback_interval) == 0
            %    fprintf("Iteration: %d, total time: %g\n", obj.iteration, cputime-obj.start_time);
                obj.callback.call(u, obj);
            end
            obj.iteration = obj.iteration + 1;
            obj.running = ~obj.termination_condition.handle(u, obj, obj.termination_condition);
        end
        function store_diff(obj, du)
           % STATE.STORE_DIFF(DU) computes the norm^2 of du and stores
           % that number in the .diffs property. This way the convergence
           % of the algorithm can be monitored
           obj.diffs = [obj.diffs, norm(du(:))^2];
        end
    end
    
end

