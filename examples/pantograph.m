function pantograph()
    close all;
    % Solves the Pantograph differential equation of the form
    % dfdt_pantograph = @(t) a*f(t) + b*f(s*t);
    % 
    
    % Parameter definition
    a = 10;  % direct derivative weight
    b = 20;  % 5 or 20 % delayed derivative weight
    s = 1/2;  % delay factor
    % Initial conditions  
    Dt_init = 1;  % Initial time period
    Dt_total = 3*60;  % 5, 150, 3*60       % Total time period runs from [0, Dt]
    % The initial function on time interval [0, Dt_init]
    f_init = @(t) 0.1 + exp(-(0.5 .* (t(:)-0.85)./0.02).^2) - 0.5*exp(-(0.5 .* (t(:)-0.80)./0.05).^2); 
    
    % Sampling
    dt = 0.01;  % Uniform sampling in time
    % Utility definitions
    nb_samples_init = round(Dt_init / dt);
    nb_samples_total = round(Dt_total / dt);
    t_range = ([1:nb_samples_total]-1) * dt;
    f_init_sampled = f_init(t_range(1:nb_samples_init).');
    
    %logMessage('Solving df(t)/dt = %0.1f f(t) + %0.1f f(%0.1ft) for f(t) on interval [%ds, %ds]', [a, b, s, Dt_init, Dt_total]);
    
    %
    % Define some functions to do continuous operations on sampled functions
    %
    function f_interpolated = interpolate(f_sampled, t)
        % Linearly interpolate the sampled function at the specified times
        % The boundary values are used when outside of the bounds
        t = t(:);
        index = floor(1 + t ./ dt);
        index = max(1, min(numel(t_range)-1, index));  % clip
        remainder = (1 + t ./ dt) - index;
        f_interpolated = (1-remainder) .* f_sampled(index) + remainder .* f_sampled(index+1);
    end
    function d = derivative(f_sampled)
        % The finite difference derivative as the average of the forward
        % and backward difference. The input is zero padded on both ends so
        % that the output vector is of the same size.
        d = diff([0; f_sampled(:); 0]);
        d = d(1:end-1) + d(2:end);
        d = d .* (0.5 / dt);
    end
    % A helper function to build a time-restricted operator from an operator
    function op_post_func = post(op)
        % A function to change a matrix multiplication operator to one that
        % only operators on the bottom right section of the matrix, i.e.
        % the part that acts on the unknown.
        function result_post = op_post(fp)
            result_post = op([zeros(nb_samples_init, 1); fp]);
            result_post = result_post(1+nb_samples_init:end);
        end
        op_post_func = @op_post;
    end
    
    % Define the linear system to solve: Hf = 0
    H = @(f_s) derivative(f_s) - a*f_s - b*interpolate(f_s, s*t_range);
    
    % Analyze using Matlab and return a function to check future tests against it
    [solution_matlab, forward_inverse_residue_func] = solve_with_matlab(H, f_init_sampled, nb_samples_total);
    
    % Rescale the problem so that ||V|| < 1
    scale = sign(-a)*abs(b/sqrt(s) / 0.9);  % |.| to make sure this is positive so that H remains accretive!
    % Define matrix multiplication operations without building the matrix in memory
    H_scaled = @(f_s) derivative(f_s)/scale - (a/scale) * f_s - (b/scale) * interpolate(f_s, s*t_range);
    G = @(f_s) f_s + (b/scale) * interpolate(f_s, s*t_range);
    LpI = @(f_s) derivative(f_s)/scale + (1 - a/scale) * f_s;  % = @(f_s) H_scaled(f_s) + G_scaled(f_s);
    V = @(f_s) f_s - G(f_s);  % G = 1-V
    % Define the scaled right-hand side
    rhs_scaled_post = -H_scaled([f_init_sampled; zeros(nb_samples_total - nb_samples_init, 1)]);
    rhs_scaled_post = rhs_scaled_post(nb_samples_init+1:end);
    
    if nb_samples_total < 1024
        % Build the matrices just to check everything is as expected
        [norm_H, max_eig_H, min_Rayleigh_Hr, max_Rayleigh_Hr, min_Rayleigh_Hi, max_Rayleigh_Hi] = check_matrix_op(post(H_scaled), nb_samples_total - nb_samples_init);
        logMessage('|H| = %0.3f, %0.3f < Re<x,Hx> < %0.3f, %0.3f < Im<x,Hx> %0.3f,  max(|eig(H)|) = %0.3f',...
            [norm_H, min_Rayleigh_Hr, max_Rayleigh_Hr, min_Rayleigh_Hi, max_Rayleigh_Hi, max_eig_H]);
        [norm_V, max_eig_V, min_Rayleigh_Vr, max_Rayleigh_Vr, min_Rayleigh_Vi, max_Rayleigh_Vi] = check_matrix_op(post(V), nb_samples_total - nb_samples_init);
        logMessage('|V| = %0.3f, %0.3f < Re<x,Vx> < %0.3f, %0.3f < Im<x,Vx> %0.3f,  max(|eig(V)|) = %0.3f',...
            [norm_V, min_Rayleigh_Vr, max_Rayleigh_Vr, min_Rayleigh_Vi, max_Rayleigh_Vi, max_eig_V]);
    end
    
    %
    % Split Richardson solver
    %
    function f_post = LpIinv_post(fp)
        % Solve the bottom-right sub-problem for H+G = L+1: post(LpI) f_post = fp
        % As LpI is a tri-diagonal Toeplitz matrix, it can be solved with a
        % linear O(N) algorithm
        below = -0.5 / (dt * scale);
        on = 1 - a/scale;
        above = -below;
        f_post = tridiagtoeplitz(below, on, above, fp);
        % Note: only invert the problem after the initial period otherwise artefacts at the boundaries!
    end
    G_post = post(G);
    prec_inv = @(x) G_post(LpIinv_post(x));
    H_post = @(x) LpI(x) - G_post(x);
    
    start_time = cputime();
    [solution_post_sr, flag, relres, iter, resvec] = splitrichardson(@LpIinv_post, G_post, rhs_scaled_post, 1e-8, 1000, zeros(nb_samples_total-nb_samples_init, 1));  %, @callback);
%     [solution_post_sr, flag, relres, iter, resvec] = bicgstab(@(x) prec_inv(H_post(x)), prec_inv(rhs_scaled_post));
    total_time = cputime() - start_time;
    %logMessage('SplitRichardson took %0.3fms', total_time*1e3);
    relative_resvec = resvec ./ norm(solution_post_sr);
    
    solution_sr = [f_init_sampled; solution_post_sr];  % prepend the initial function
    
    %logMessage('SplitRichardson forward residue %0.6f and inverse residue %0.6f.', forward_inverse_residue_func(solution_sr));
    
    % Double check M
    M = @(f_s) G_post(LpIinv_post(G_post(f_s))) + f_s - G_post(f_s);
    if nb_samples_total < 1024
        % Build the matrices just to check everything is as expected
        [norm_M, max_eig_M, min_Rayleigh_Mr, max_Rayleigh_Mr, min_Rayleigh_Mi, max_Rayleigh_Mi] = check_matrix_op(M, nb_samples_total - nb_samples_init);
        %logMessage('|M| = %0.3f, %0.3f Re<x,Mx> < %0.3f, %0.3f Im<x,Mx> %0.3f,  max(|eig(M)|) = %0.3f',...
        %    [norm_M, min_Rayleigh_Mr, max_Rayleigh_Mr, min_Rayleigh_Mi, max_Rayleigh_Mi, max_eig_M]);
    end
    
    %
    % Display results
    %
    fig = figure('Position', [50, 50, 1024, 1024], 'Color', [1, 1, 1]);
    ax(1) = subplot(2, 1, 1);
    plot(t_range, solution_sr*0, 'Color', [0, 0, 0], 'LineWidth', 0.5); hold on;
    plot([1, 1], [min(solution_sr), max(solution_sr)].*1.05, 'Color', [0, 0, 0], 'LineWidth', 0.5); hold on;
    plot(t_range, solution_matlab, 'Color', [0.8, 0, 0], 'LineWidth', 3); hold on;
    plot(t_range, solution_sr, 'Color', [0, 0.75, 0], 'LineWidth', 1.5); hold on;
%     plot(t_range(1:nb_samples_init), solution_matlab(1:nb_samples_init), 'Color', [0.8, 0, 0], 'LineWidth', 6);
%     plot(t_range(1:nb_samples_init), solution_sr(1:nb_samples_init), 'Color', [0, 0.75, 0], 'LineWidth', 3);
    box off;
    xlim(t_range([1, end]));
    ylim([min(solution_sr), max(solution_sr)].*1.05);
    xlabel('t [s]', 'FontSize', 24);
    ylabel('f(t)  [a.u.]', 'FontSize', 24);
    legend({'', '', 'LU inverse solver', 'Split-Richardson'});
    
    ax(2) = subplot(2, 1, 2);
    semilogy(relative_resvec, 'LineWidth', 3);
    xlim([1, numel(relative_resvec)]);
    ylim([1e-8, 1]);
    box off;
    xlabel('iteration #', 'FontSize', 24)
    ylabel('relative residue Split Richardson', 'FontSize', 24)
    set(ax, 'FontSize', 18, 'LineWidth', 2);
    
    plot2svg(sprintf('pantograph%d.svg', b), fig)
end

function [norm_H, max_eig_H, min_Rayleigh_Hr, max_Rayleigh_Hr, min_Rayleigh_Hi, max_Rayleigh_Hi] = check_matrix_op(H, nb_samples_total)
    % Build the matrices column-by-column just to check everything is correct
    H_sampled = zeros(nb_samples_total, nb_samples_total);
    tmp = zeros(nb_samples_total, 1);
    for col_idx = [1:nb_samples_total]
        tmp(col_idx) = 1;
        H_sampled(:, col_idx) = H(tmp);
        tmp(col_idx) = 0;
    end
    clear tmp

    % Analyze the matrix
    norm_H = norm(H_sampled);
    max_eig_H = max(abs(eig(H_sampled)));
    Dr = real(eig(H_sampled + H_sampled')) ./ 2;
    min_Rayleigh_Hr = min(Dr);
    max_Rayleigh_Hr = max(Dr);
    Di = real(eig((1i*H_sampled) + (1i*H_sampled)')) ./ 2;
    min_Rayleigh_Hi = min(Di);
    max_Rayleigh_Hi = max(Di);
end

function [solution_matlab, residue_function] = solve_with_matlab(H, f_init_sampled, nb_samples_total)
    nb_samples_init = numel(f_init_sampled);
    H_sampled = zeros(nb_samples_total, nb_samples_total, 'single');
    tmp = zeros(nb_samples_total, 1);
    for col_idx = [1:nb_samples_total]
        tmp(col_idx) = 1;
        H_sampled(:, col_idx) = H(tmp);
        tmp(col_idx) = 0;
    end
    clear tmp
    
    H_sampled_post_pre = H_sampled(nb_samples_init+1:end, 1:nb_samples_init);
    H_sampled_post_post = H_sampled(nb_samples_init+1:end, nb_samples_init+1:end);
    % Define the right-hand side
    rhs_post = -H([f_init_sampled; zeros(nb_samples_total - nb_samples_init, 1)]);
    rhs_post = rhs_post(nb_samples_init+1:end);
    
    start_time = cputime();
    solution_matlab_post = H_sampled_post_post \ rhs_post;
    total_time = cputime() - start_time;
    %logMessage('MATLAB took %0.3fms', total_time*1e3);
    solution_matlab = [f_init_sampled; solution_matlab_post];  % prepend the initial function
    
    % Define a residue function to check against later
    function forward_inverse_residue = calc_residues(solution)
        forward_res = norm([H_sampled_post_pre, H_sampled_post_post] * solution) / norm(solution_matlab);
        inverse_res = norm(solution - solution_matlab) / norm(solution_matlab);
        forward_inverse_residue = [forward_res, inverse_res];
    end
    residue_function = @calc_residues;
    
    %logMessage('MATLAB forward residue %0.6f and inverse residue %0.6f', residue_function(solution_matlab));    
end

function x = tridiagtoeplitz(below, on, above, y)
    % Solves Ax = y where A is a Toeplitz tri-diagonal matrix of the form:
    % [ on   above   0     0     0  ]
    % [below  on   above   0     0  ]
    % [ 0    below  on   above   0  ]
    % [ 0      0 below     on  above]
    % [ 0      0    0    below   on ]
    % https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    %
    y = y(:);
    n = size(y, 1);
    % Forward sweep
    c_prime = zeros(n, 1);
    d_prime = zeros(n, 1);
    c_prev = 0.0;
    d_prev = 0.0;
    for idx = [1:n]
        d_prev = (y(idx) - below * d_prev) / (on - below * c_prev);
        d_prime(idx) = d_prev;
        c_prev = above / (on - below * c_prev);
        c_prime(idx) = c_prev;
    end
    clear c_prev d_prev
    
    % Back substitution
    x = d_prime;
    for idx = [n-1:-1:1]
        x(idx) = x(idx) - c_prime(idx) * x(idx+1);
    end
end