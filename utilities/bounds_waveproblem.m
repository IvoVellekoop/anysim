% computes the combination of A and V for which |M| is maximized
% assumes that L and V are both accretive
% conclusion: max |M| is achieved when |imag V| is maximized (largest
% refractive index perturbation) and real A = 0, imag A = Amin
% This point corresponds to L imaginary.
% To encounter this in practice, we would have a homogeneous medium
% filled at a constant refractive index (V=Vmax almost everywhere)

clear; close all;
opt.Vmax = 0.95;
opt.Amin = 1E-1;
opt.Amax = 1E1;
opt.N = 512;

alpha_r = 0:0.1:1;
Nalpha = length(alpha_r);
Mmax = zeros(Nalpha, 1);
for a_i = 1:Nalpha
    Mmax(a_i) = plot_all(alpha_r(a_i), opt);
    xlim([.95 1])
    ylim([-1E-1 1E-1])
    viscircles([0,0], 1)
    pause(1);
end
figure; plot(alpha_r, 1-abs(Mmax));

function Mval = plot_all(alpha, opt)
    % A spans half of complex plane, with small semicircle cut out around Amin
    A1 = linspace(1.0i * opt.Amin, 1.0i * opt.Amax, opt.N).';
    A2 = opt.Amax * exp(1i * linspace(pi/2, -pi/2, opt.N)).';
    A3 = linspace(-1.0i * opt.Amax, -1.0i * opt.Amin, opt.N).';
    A4 = opt.Amin * exp(1i * linspace(-pi/2, pi/2, opt.N)).';
    A = [A1; A2; A3; A4];
    
    V1 = opt.Vmax * linspace(-1.0i, 1.0i, opt.N);
    V2 = opt.Vmax * exp(1i * linspace(pi/2, -pi/2, opt.N));
    V = [V1, V2];
    
    M = calc_M(alpha, A, V);
    
    plot(M(:));
    [Mmax, ind] = max(abs(M(:)));
    Mval = M(ind);
    [a_i, v_i] = ind2sub(size(M), ind);
    hold on;
    plot(M(ind), '+');
    hold off;
    fprintf('min 1-|M| %g at A=%g + %f i, V=%f + i%f i\n', 1-Mmax, A(a_i), imag(A(a_i)), real(V(v_i)), imag(V(v_i)));
end

function M = calc_M(alpha, A, V)
    M = 1 - alpha * (1-V).*(1-1./(A+1-V).*(1-V));
end