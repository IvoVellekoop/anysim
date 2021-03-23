function terminate = tc_relative_error(u, state, opt)
%TC_RELATIVE_ERROR Simple termination condition
    reldiff = state.diffs / max(state.diffs);
    terminate = reldiff < opt.relative_error_limit;
end

