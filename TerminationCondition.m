classdef TerminationCondition
    %TERMINATIONCONDITION Default termination condition
    %   Terminates the algorithm if any of the following conditions is true:
    %       * The maximum number of iterations (specified in
    %         OPT.iteration_count) is reached
    %       * The 'magnitude' of the change in solution 'u'
    %         is below a certain absolute limit (specified in 
    %         OPT.absolute_error). The strength of the update
    %         is calculated by the function OPT.quantify_update.
    %         For grid-based simulations, by default this function
    %         crops 'u' to the roi and computes the norm squared.
    %       * The 'magnitude' of the change in solution 'u'
    %         is below a certain relative limit (specified in 
    %         OPT.relative_error). The value is relative to the
    %         largest update magnitude encountere so far. (should be the
    %         1st!)
    %
    properties
        iteration_count % maximum interation count
        absolute_limit %
        relative_limit %
    end
    
    methods
        function obj = TerminationCondition(sim, opt)
            %TERMINATIONCONDITION Don't call this function directly
            % see sample scripts.
            default.iteration_count = 1E4;
            default.absolute_limit = 1E-12;
            default.relative_limit = 1E-6;
            default.show_boundaries = false;
            default.show_convergence = false;
            opt = set_defaults(default, opt);
            obj.iteration_count = opt.iteration_count;
            obj.absolute_limit = opt.absolute_limit;
            obj.relative_limit = opt.relative_limit;
        end
        
        function terminate = call(obj, state)
            current_diff = state.diffs(end);
            terminate = state.iteration >= obj.iteration_count ||...
                current_diff <= obj.absolute_limit || ...
                current_diff / max(state.diffs) <= obj.relative_limit;
        end
    end
end

