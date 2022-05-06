function inspect_sim(sim)
    arguments
        sim GridSim
    end
%INSPECT_SIM(SIM) Verifies an AnySim simulation object
%   Checks if A is accretive
%   Checks if ||V|| < 1
%   Tests if A is Hermitian or skew-Hermitian or neither
%
    N = sim.grid.N_u;
    propagator = full_matrix(sim.propagator, N);
    L1 = inv(propagator); % L+1
    B = full_matrix(sim.medium, N);
    A = L1 - B;

    % verify that norm V = 1-B < 1
    V = B - eye(size(B, 1));
    nV = norm(V);
    fprintf("||V|| = %f\n", nV);
    if nV > 1
        warning('||V|| > 1');
    elseif nV > sim.opt.V_max
        warning('||V|| > V_max');
    end

    % verify that A is accretive
    eA = eig(A + A');
    if any(real(eA) < 0)
        warning(['A is not accretive, Re A = ' num2str(min(real(eA)))]);
    end

    % plot eigenvalues of A
    figure;
    title('Eigenvalues');
    eA = eig(A);
    eL = eig(L1);
    plot(real(eA), imag(eA), '+', real(eL), imag(eL), '+');
    legend('A', 'L');
    niA = norm(inv(A));
    fprintf("||A^{-1}|| = %f\n", niA);
    fprintf("||A^{-1}||||V|| = %f\n", niA * nV);

    test_hermitian(V, "V");
    test_hermitian(A, "A");
    test_hermitian(L1 - eye(size(L1, 1)), "L");
end

function test_hermitian(M, name)
    reM = 0.5 * (M + M');
    imM = 0.5 * (M - M');
    nM = norm(M(:));
    symRel = norm(reM(:)) / nM;
    asymRel = norm(imM(:)) / nM;

    if symRel < 1E-2
        fprintf("%s is skew-Hermitian to within %f\n", name, symRel);
    elseif asymRel < 1E-3
        fprintf("%s is (approximately) Hermitian to within %f\n", name, asymRel);
    else
        fprintf("%s is neither Hermitian nor skew-Hermitian\n", name);
    end
end

