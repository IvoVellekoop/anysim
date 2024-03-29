classdef State < dynamicprops
    %STATE Holds state variables (variables that change while running the
    %   algorithm.)
    %
    %   (c) 2019. Ivo Vellekoop    
    properties (SetAccess = private)
        iteration;
        start_time;
        run_time;
        termination_condition_interval;
        termination_condition;
        callback_interval;
        callback;
        running;
        residuals;
        residual_its; % iteration numbers at which residuals were reported
        normb; % norm of the preconditioned source (for normalizing 'residual')
    end
    properties
        source = []; % for use in matlab-style iterative schemes 
        internal_iteration; % number of internal iterations for two-step preconditioners
        internal_iteration_failed = false; % flag to indicate whether internal iteration failed to converge
    end
    
    methods
        function obj = State(sim, opt)
            obj.callback_interval = opt.callback_interval;
            obj.callback = opt.callback.prepare(sim);
            obj.termination_condition_interval = opt.termination_condition_interval; % how often to analyze the update du (for termination condition and/or visual feedback)
            obj.termination_condition = opt.termination_condition.prepare(sim);
            obj.internal_iteration = 0;
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
                nr = norm(r(:));
                obj.residuals = [obj.residuals, nr/obj.normb];
                obj.residual_its = [obj.residual_its, obj.iteration];
                obj.running = ~obj.termination_condition.call(obj);
                if ~isfinite(nr)
                    if obj.running
                        %warning("Terminating simulation, found non-finite residual %f", nr);
                        error("Terminating simulation, found non-finite residual %f", nr);
                        obj.running = false;
                    end
                end
            end
            if mod(obj.iteration-1, obj.callback_interval) == 0 && ~isempty(obj.callback)
                obj.callback.call(u, r, obj);
            end
            obj.iteration = obj.iteration + 1;
        end
        function reset(obj)
            obj.iteration = 1; %current iteration number
            obj.start_time = tic; % using 'tic' and not 'cputime' because we need to include the gputime
            obj.running = true;
            obj.normb = [];
            obj.residuals = [];
            obj.residual_its = [];
        end
        function finalize(obj)
            obj.run_time = toc(obj.start_time);
        end
    end
    
end

