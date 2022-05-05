function [L, GL] = simulation_eigenvalues(sim)
%% Computes eigenvalues of the preconditionaed and non-predonditioned operator
%% Note: requires the option sim.opt.forward_operator = true
    N = prod(sim.grid.N_u);
    
    GL = eig(full_matrix(sim.preconditioned, N));
    if sim.opt.forward_operator
        L = eig(full_matrix(sim.operator, N));
    else
        L = [];
    end
end

