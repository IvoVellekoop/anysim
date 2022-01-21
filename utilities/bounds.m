clear; close all;
% L and V both accretive
Vmax = 0.95;
alpha = 1;
Lmax = 1E2;
NL = 1025;
NV = 128;

L1 = Lmax * linspace(-1.0i, 1.0i, NL).';
L2 = Lmax * exp(1i * linspace(pi/2, -pi/2, NV)).';
L = [L1; L2];

V1 = Vmax * linspace(-1.0i, 1.0i, NV);
V2 = Vmax * exp(1i * linspace(pi/2, -pi/2, NV));
V = [V1, V2];

M = calc_M(alpha, L, V);

plot(M(:));
[Mmax, ind] = max(abs(M(:)));
[l_i, v_i] = ind2sub(size(M), ind);
hold on;
plot(M(ind), '+');
hold off;
fprintf('max M %f at L=%f, V=%f', Mmax, L(l_i), V(v_i));

function M = calc_M(alpha, L, V)
    Gamma_inv = alpha * (1-V) ./ (L+1); 
    M = 1 - Gamma_inv .* (L+V);
end