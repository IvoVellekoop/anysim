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
        iteration_count (1,1) double = 1E4
        absolute_limit (1,1) double = 1E-12
        relative_limit (1,1) double = 1E-3
    end
    
    methods
        function obj = TerminationCondition(opt)
            arguments
                opt.?TerminationCondition
            end
            obj = copy_properties(obj, opt);
        end

        function obj = prepare(obj, sim)
        end
        
        function terminate = call(obj, state)
            current_residual = state.residuals(end);
            terminate = state.iteration >= obj.iteration_count ||...
                current_residual <= obj.absolute_limit || ...
                current_residual / max(state.residuals) <= obj.relative_limit;
        end
    end
end

