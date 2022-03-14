function [L, GL] = simulation_eigenvalues(sim)
%% Computes eigenvalues of the preconditionaed and non-predonditioned operator
%% Note: requires the option sim.opt.forward_operator = true
N = prod(sim.grid.N);

GL = eig(full_matrix(sim.preconditioned));
if sim.opt.forward_operator
    L = eig(full_matrix(sim.operator));
else
    L = [];
end

    function M = full_matrix(A)
        M = zeros(N, N);
        b = zeros(N, 1);
        b(1) = 1;
        for n=1:N
            M(:,n) = A(b);
            b = circshift(b, [1,0]);
        end
    end
end

