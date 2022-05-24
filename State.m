classdef State < dynamicprops
    %STATE Holds state variables (variables that change while running the
    %   algorithm.)
    %
    %   (c) 2019. Ivo Vellekoop    
    properties (SetAccess = private)
        iteration;
        start_time;
        end_time;
        run_time;
        termination_condition_interval;
        termination_condition;
        callback_interval;
        callback;
        running;
        residuals;
        residual_its; % iteration numbers at which residuals were reported
        normb; % norm of the source (for normalizing 'residual')
    end
    properties
        source = []; % for use in matlab-style iterative schemes 
    end
    
    methods
        function obj = State(sim, opt)
            obj.callback_interval = opt.callback_interval;
            obj.callback = opt.callback.prepare(sim);
            obj.termination_condition_interval = opt.termination_condition_interval; % how often to analyze the update du (for termination condition and/or visual feedback)
            obj.termination_condition = opt.termination_condition.prepare(sim);
            obj.reset();
        end
        function next(obj, u, r)
            % when running AnySim, r is the residual b-Ax
            % when running a MATLAB iterative algorithm, r equals Ax only
            % to compute the residual, we need to use the stored source
            % value
            if ~isempty(obj.source) && (...
                    mod(obj.iteration-1, obj.termination_condition_interval) == 0 ||...
                    (mod(obj.iteration-1, obj.callback_interval) == 0 && ~isempty(obj.callback)))
                r = reshape(obj.source, size(r)) - r;
            end
            
            if mod(obj.iteration-1, obj.termination_condition_interval) == 0
                if (obj.iteration == 1)
                    obj.normb = norm(r(:));
                end
                obj.residuals = [obj.residuals, norm(r(:))/obj.normb];
                obj.residual_its = [obj.residual_its, obj.iteration];
                obj.running = ~obj.termination_condition.call(obj) && isfinite(obj.residuals(end));
            end
            if mod(obj.iteration-1, obj.callback_interval) == 0 && ~isempty(obj.callback)
                obj.callback.call(u, r, obj);
            end
            obj.iteration = obj.iteration + 1;
        end
        function reset(obj)
            obj.iteration = 1; %current iteration number
            obj.start_time = cputime; %only measure time actually used by MATLAB
            obj.running = true;
            obj.normb = [];
            obj.residuals = [];
            obj.residual_its = [];
        end
        function finalize(obj)
            obj.end_time = cputime;
            obj.run_time = obj.end_time - obj.start_time;
        end
    end
    
end

