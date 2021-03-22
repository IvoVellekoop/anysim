function terminate = tc_fixed_iteration_count(u, state, opt) %#ok<INUSL>
%TC_FIXED_ITERATION_COUNT Simple termination condition
    terminate = state.iteration > opt.iteration_count;
end

