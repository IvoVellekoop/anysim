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
    M = eye(size(propagator, 1)) - full_matrix(sim.preconditioned, N);
    L1 = inv(propagator); % L+1
    B = full_matrix(sim.medium, N);
    A = L1 - B;
    L = L1 - eye(size(L1, 1));
    A11 = A(1:2:end, 1:2:end);
    A12 = A(1:2:end, 2:2:end);
    A21 = A(2:2:end, 1:2:end);
    A22 = A(2:2:end, 2:2:end);

    % verify that norm V = 1-B < 1
    V = B - eye(size(B, 1));
    nV = norm(V);
    fprintf("||V|| = %f\n", nV);
    if nV > 1
        warning('||V|| > 1');
    elseif nV > sim.opt.V_max
        warning('||V|| > V_max');
    end

    nM = norm(M);
    if nM > 1
        eM = eig(M);
        warning('||M|| = %f > 1. Note: max |λ_M| = %f', nM, max(abs(eM)));
    end
    % compute eigenvalues of the operators
    eAh = eig(0.5 * (A + A'));
    eA = eig(A);
    eL = eig(L);

    % verify that A is accretive
    if any(real(eAh) < 0)
        if any(real(eA) < 0)
            warning('A has negative eigenvalues, min λ_A = %f', min(real(eA)));
        else
            warning("A is not accretive, but all eigenvalues are positive, Re A = %f", min(real(eAh)));
        end
    end

    % plot eigenvalues of A
    figure;
    title('Eigenvalues');
    plot(real(eA), imag(eA), '+', real(eL), imag(eL), '+');
    legend('A', 'L');
    niA = norm(inv(A));
    fprintf("||A^{-1}|| = %f\n", niA);
    fprintf("||A^{-1}||||V|| = %f\n", niA * nV);

    test_hermitian(V, "V");
    test_hermitian(A, "A");
    test_hermitian(L, "L");
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

