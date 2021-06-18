function compareSplitRichardson(nb_vars)
    close all;
    if nargin < 1
        nb_vars = 2^16;
    end
    
    rng = RandStream('mt19937ar', 'Seed', 1);
    
    nb_vectors = min(nb_vars, 10);
    nb_trials = 100;
    maxit = 100;
    tol = eps('single');  % Smaller than achievable
    restarts = [];  %min(10, nb_vars);  %default for gmres
    
    function result = rand_vec()
        result = (rng.randn(nb_vars, 1) + 1i * rng.randn(nb_vars, 1)) ./ sqrt(2);
    end

    evals = [0, 0, 0];
    function result = H_func_counter(x, L_diag, V_diag, V_max)
        evals(1) = evals(1) + 1;
        result = ifft(fft(x) .* L_diag) + (V_diag + V_max) .* x;
    end

    function result = HpGinv_func_counter(x, L_diag, V_max)
        evals(2) = evals(2) + 1;
        result = ifft(fft(x) ./ (L_diag + (1 + V_max)));
    end

    function result = G_func_counter(x, V_diag)
        evals(3) = evals(3) + 1;
        result = x - V_diag .* x;
    end

    function values = transp(values, transp_flag)
        if strcmpi(transp_flag, 'transp')
            values = conj(values);
        end
    end

    function cont = generic_callback(iter, relres, x, dx, H_diag, V_diag, V_max, b, tol_x_b)
        cont = mod(iter+1, 1) > 0 || norm(H_func_counter(x, H_diag, V_diag, V_max) - b) >= tol_x_b;
    end

    % Auxiliary function
    function memused = show_memory(label)
%         [userview, systemview] = memory();
%         memused = userview.MemUsedMATLAB;
%         logMessage([label, ' %0.9f GiB'], memused / 2^30);
    end
    
    max_nb_algorithms = 7;
    
    psi_err = NaN(nb_trials, nb_vectors, max_nb_algorithms);
    psi_rhs_err = NaN(nb_trials, nb_vectors, max_nb_algorithms);
    times = NaN(nb_trials, nb_vectors, max_nb_algorithms);
    evaluations = NaN(nb_trials, nb_vectors, max_nb_algorithms, numel(evals));
    legends = {};
    for trial_idx = 1:nb_trials
        logMessage('Trial %d/%d for %dx%d with %d random vectors...', [trial_idx, nb_trials, nb_vars, nb_vars, nb_vectors]);
        % Pick random matrices
        V_diag = rand_vec();
        V_max = rng.rand();
        V_diag = V_diag .* (V_max / max(abs(V_diag)));
%         V_diag = abs(real(V_diag)) + 1i * imag(V_diag);  %TODO: remove
        L_diag = 10 * rand_vec();
        L_diag = abs(real(L_diag)) + 1i * imag(L_diag);
%         L_diag = 5 * (([1:numel(V_diag)].' - 1) ./ numel(V_diag)).^2;  % TODO: remove
        % Scale the problem
        scale = 0.5 ./ max(abs(V_diag));
        L_diag = single(L_diag .* scale);
        V_diag = single(V_diag .* scale);
        
        % Define functions
        H_func = @(x) H_func_counter(x, L_diag, V_diag, V_max);
        G_func = @(x) G_func_counter(x, V_diag);
        % Define the preconditioner
        % HpG_func = @(x) ifft(fft(x) .* (H_diag)) + (1 + V_max) .* x;
        % HpG_func = @(x) ifft(fft(x) .* (H_diag + 1 + V_max));
        HpGinv_func = @(x) HpGinv_func_counter(x, L_diag, V_max);
        prec_inv_func = @(x) G_func(HpGinv_func(x));
        
        % Some algorithms need the transposes as well:
        H_func_with_transpose = @(x, t) H_func_counter(x, transp(L_diag, t), transp(V_diag, t), V_max);
        prec_inv_func_with_transpose = @(x, t) HpGinv_func_counter(G_func_counter(x, transp(V_diag, t)), transp(L_diag, t), V_max);
        
%         % TEST
%         H = zeros(nb_vars, nb_vars);
%         H2 = H;
%         Ht = H;
%         for col_idx = [1:nb_vars]
%             x = zeros(nb_vars, 1);
%             x(col_idx) = 1;
%             H(:, col_idx) = H_func(x);
%             H2(:, col_idx) = H_func_with_transpose(x, 'notransp');
%             Ht(:, col_idx) = H_func_with_transpose(x, 'transp');
%         end

        for vector_idx = 1:nb_vectors
            logMessage('Trial %d/%d for %dx%d with random vector %d/%d...', [trial_idx, nb_trials, nb_vars, nb_vars, vector_idx, nb_vectors]);
            % Define a problem to solve
            ground_truth = rand_vec();
            b = H_func(ground_truth);
            norm_b = norm(b);

            %
            % solve it
            %
            plot_idx = 0;
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             psi = H \ b;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H*psi - b) / norm_b;
%             leg{plot_idx} = 'Matlab \\';
% 
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             psi = (prec_inv * H) \ (prec_inv * b);
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H*psi - b) / norm_b;
%             leg{plot_idx} = 'Matlab \\+\Gamma';

            plot_idx = plot_idx + 1;
            evals = [0, 0, 0];
            show_memory('before splitrichardson');
            tic();
%             [psi, flag, relres, iter, resvec] = splitrichardson(HpGinv_func, G_func, b, tol, maxit, [], ...
%                 @(iter, relres, x, dx) generic_callback(iter, relres, x, dx, L_diag, V_diag, V_max, b, tol*norm_b));
            [psi, flag, relres, iter, resvec] = splitrichardson(HpGinv_func, G_func, b, 0, maxit);  % Note the x4 !
            evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
            times(trial_idx, vector_idx, plot_idx) = toc();
            psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
            psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
            legends{plot_idx} = 'splitrichardson';
            show_memory('after splitrichardson');
            
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             [psi, flag, relres, iter] = gmres(H_func, b, restarts, tol, min(maxit, nb_vars));
%             evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
%             leg{plot_idx} = 'GMRES';

%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             [psi, flag, relres, iter, resvec] = gmres(H_func, b, restarts, tol, min([maxit, nb_vars, floor(1e8 / nb_vars)]), prec_inv_func);
%             evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
%             leg{plot_idx} = 'GMRES+\Gamma';
% 
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             [psi, flag, relres, iter] = bicg(H_func_with_transpose, b, tol, min(maxit, nb_vars), prec_inv_func_with_transpose);
%             evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx,plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
%             leg{plot_idx} = 'bicg+\Gamma';
% 
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             [psi, flag, relres, iter] = bicgstab(H_func, b, tol, min(maxit, nb_vars));
%             evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx, plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
%             leg{plot_idx} = 'bicgstab';

            plot_idx = plot_idx + 1;
            evals = [0, 0, 0];
            show_memory('before bicgstab+Gamma');
            tic();
            warning('off', 'MATLAB:bicgstab:tooSmallTolerance')
            [psi, flag, relres, iter] = bicgstab(H_func, b, tol, min(maxit, nb_vars), prec_inv_func);
            evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
            times(trial_idx, vector_idx, plot_idx) = toc();
            psi_err(trial_idx, vector_idx,plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
            psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
            legends{plot_idx} = 'bicgstab+\Gamma';
            show_memory('after bicgstab+Gamma');
            
%             plot_idx = plot_idx + 1;
%             evals = [0, 0, 0];
%             tic();
%             [psi, flag, relres, iter] = bicgstabl(H_func, b, tol, min(maxit, nb_vars), prec_inv_func);
%             evaluations(trial_idx, vector_idx, plot_idx, :) = evals;
%             times(trial_idx, vector_idx, plot_idx) = toc();
%             psi_err(trial_idx, vector_idx,plot_idx) = norm(psi - ground_truth) / norm(ground_truth);
%             psi_rhs_err(trial_idx, vector_idx, plot_idx) = norm(H_func(psi) - b) / norm_b;
%             leg{plot_idx} = 'bicgstab(l)+\Gamma';
        end
    end
    
    % Calculate statistics
    times = reshape(times, nb_trials * nb_vectors, []);
    psi_err = reshape(psi_err, nb_trials * nb_vectors, []);
    psi_rhs_err = reshape(psi_rhs_err, nb_trials * nb_vectors, []);
    evaluations = reshape(evaluations, nb_trials * nb_vectors, [], numel(evals));
    
    times_median = nanmedian(times);
    times_std = nanstd(times);
    times_max = nanmax(times);
    psi_err_median = nanmedian(psi_err);
    psi_err_std = nanstd(psi_err);
    psi_err_max = nanmax(psi_err);
    psi_rhs_err_median = nanmedian(psi_rhs_err);
    psi_rhs_err_std = nanstd(psi_rhs_err);
    psi_rhs_err_max = nanmax(psi_rhs_err);
    evaluations_median = nanmedian(evaluations, 1);
    evaluations_median = reshape(evaluations_median, [], numel(evals)).';
    evaluations_std = nanstd(evaluations, 1);
    evaluations_std = reshape(evaluations_std, [], numel(evals)).';
    evaluations_max = nanmax(evaluations, [], 1);
    evaluations_max = reshape(evaluations_max, [], numel(evals)).';
    
    logMessage(['median    err/tol: ', sprintf('%0.9f, ', psi_rhs_err_median/tol)]);
    logMessage(['median conv evals: ', sprintf('%0.0f, ', sum(evaluations_median(1:2, :), 1))]);
    logMessage(['median       time: ', sprintf('%0.9f, ', times_median)]);
    logMessage(['   max    err/tol: ', sprintf('%0.9f, ', psi_rhs_err_max/tol)]);
    logMessage(['   max conv evals: ', sprintf('%0.0f, ', sum(evaluations_max(1:2, :), 1))]);
    logMessage(['   max       time: ', sprintf('%0.9f, ', times_max)]);
    
    % Display 
    colors = [0, 0.5, 0; 1, 0, 0; 0, 0, 1; 0, 0.5, 1.0; 1.0, 0.5, 0.0; 1.0, 0.0, 0.5; 0.5, 1.0, 0.5; 0.5, 0.5, 0.5];
    line_widths = [3, 2 * ones(1, size(colors, 1)-1)];
    function plot_hist(values, integer)
        if nargin < 2 || isempty(integer)
            integer = false;
        end
        if ~integer
            nb_bins = min(max(10, size(values,1) / 10), 50);
            min_value = quantile(values(~isnan(values)), 0.001);
            max_value = quantile(values(~isnan(values)), 0.999);
            edges = ([0:nb_bins] ./ nb_bins) .* (max_value - min_value) + min_value;
            edges(1) = min(values(:));  % Make the first bin sufficiently large to include all
            edges(end) = max(values(:));  % Make the last bin sufficiently large to include all
        else
            edges = [-0.5:max(values(:))+0.5];
            nb_bins = numel(edges) - 1;
        end
        nb_hists = size(values, 2);
        counts = zeros([nb_bins, nb_hists]);
        for hist_idx = [1:nb_hists]
            [counts(:, hist_idx), edges] = histcounts(values(:, hist_idx), edges);
            centers = (edges(1:end-1) + edges(2:end)) ./ 2;
        end
        hold off;
        for hist_idx = [1:nb_hists]
            plot(edges(floor(1.5:0.5:end)), counts(floor(1:0.5:end+0.5), hist_idx), 'Color', colors(hist_idx, :), 'LineWidth', line_widths(hist_idx)); hold on;
        end
        for hist_idx = [1:nb_hists]
            scatter(centers, counts(:, hist_idx), 25, colors(hist_idx, :), 'filled'); hold on;
        end
        hold off;
    end
    
    fig = figure('Position', [50, 50, 1024, 768], 'Name', sprintf('%dx%d', [nb_vars, nb_vars]));
    axs(1) = subplot(2, 4, 1);
    plot_hist(psi_err);
    xlabel('\epsilon_\psi'); ylabel('#');
    axs(1) = subplot(2, 4, 2);
    plot_hist(psi_rhs_err);
    xlabel('\epsilon_b'); ylabel('#');
    axs(3) = subplot(2, 4, 3);
    plot_hist(times);
    xlabel('t  [s]'); ylabel('#');
%     axs(4) = subplot(2, 4, 4);
%     plot_hist(evaluations(:,:,1) + 0*cat(2, zeros(size(evaluations, 1), 1, 1), evaluations(:,2:end,2)), true);
%     xlabel('convs - criterion checks'); ylabel('#');
    axs(5) = subplot(2, 4, 5);
    plot_hist(evaluations(:,:,1), true);
    xlabel('H evals'); ylabel('#');
    axs(6) = subplot(2, 4, 6);
    plot_hist(evaluations(:,:,2), true);
    xlabel('inv(H+G) evals'); ylabel('#');
    axs(7) = subplot(2, 4, 7);
    plot_hist(evaluations(:,:,3), true);
    xlabel('G evals'); ylabel('#');
    axs(8) = subplot(2, 4, 8);
    plot_hist(evaluations(:,:,1) + evaluations(:,:,2), true);
    xlabel('convs = H + inv(H+G)'); ylabel('#');
    legend(legends);
    
    drawnow();
    savefig(fig, 'compareSplitRichardson_result.fig')
 end