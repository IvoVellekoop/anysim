function terminate = tc_relative_error(u, state, opt)
%TC_RELATIVE_ERROR Simple termination condition
    if ~isprop(state, 'absolute_error_limit')
        state.addprop('absolute_error_limit');
        state.absolute_error_limit = max(abs(u(:))) * opt.relative_error_limit;
    end
    if ~isprop(state, 'stored_u')
        state.addprop('stored_u');
        terminate = false;
    else
        absolute_error = max(abs(u(:) - state.stored_u(:)));
        terminate = absolute_error < state.absolute_error_limit;
    end
    state.stored_u = u;
end

