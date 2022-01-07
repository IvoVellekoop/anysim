function [L, GL] = simulation_eigenvalues(sim)
%% Computes eigenvalues of the preconditionaed and non-predonditioned operator
%% Note: requires the option sim.opt.forward_operator = true

N = prod(sim.grid.N);
GA = zeros(N(1), N(1));
A = zeros(N(1), N(1));
y = zeros(N(1), 1);
y(1) = 1;
for n=1:N(1)
    GA(:,n) = sim.preconditioned(y);
    A(:,n) = sim.operator(y);
    y = circshift(y, [1,0]);
end

%%
GL = eig(GA);
L = eig(A);
%fprintf('||A^{-1}|| = %f\n', abs(1./min(E)))
%fprintf('inf ||A^{-1}|| = %f\n', abs(1./max(E)))