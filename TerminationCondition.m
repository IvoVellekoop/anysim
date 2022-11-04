classdef TerminationCondition
    %TERMINATIONCONDITION Default termination condition
    %   Terminates the algorithm if any of the following conditions is true:
    %       * The maximum number of iterations (specified in
    %         OPT.iteration_count) is reached
    %       * The residual is below a certain absolute limit (specified in 
    %         OPT.absolute_error). The strength of the update
    %         is calculated by the function OPT.quantify_update.
    %         For grid-based simulations, by default this function
    %         crops 'u' to the roi and computes the norm squared.
    %       * The residual is below a certain relative limit (specified in 
    %         OPT.tolerance). The value is relative to the norm of
    %         the source vector. This is the typical 'tolerance' used in
    %         MATLAB builtin functions.
    %       * The residual is above a certain relative limit, or it is a NaN.
    %         this may happen if the system to be solved is not accretive,
    %         or if there is a bug in the implementation.
    %
    properties
        iteration_count (1,1) double = 1E4
        absolute_limit (1,1) double = 1E-12
        tolerance (1,1) double = 1E-3
        divergence_limit (1,1) double = 1E12
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
                current_residual / state.normb <= obj.tolerance || ...
                current_residual / state.normb >= obj.divergence_limit;
        end
    end
end

