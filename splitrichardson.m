function [x, flag, relres, iter, resvec] = splitrichardson(LpIinv, ImV, b, tol, maxit, x0, callback, display_progress)
  %
  % [x, flag, relres, iter, resvec] = splitrichardson(LpIinv, ImV, b, tol, maxit, x0, callback, display_progress)
  %
  % Solves the (non-symmetric) linear system Hx = (L+V)x = b iteratively
  % for x when the real part of H=L+V is positive definite, i.e. real(x'*H*x) > 0;
  % and V  is a contraction, i.e. norm(V) < 1. The latter can be ensured 
  % for any bound linear system by scaling both sides appropriately.
  % Similarly, when the Rayleigh quotient of (L+V) is contained within any
  % half of the complex plane, multiplication by exp(-1i*phi) for some
  % angle phi can ensure that the real part of its Rayleigh quotient is
  % positive.
  %
  % This method is computationally efficient when:
  % - The inverse, (L+I)\y, can be calculated efficiently, where I is the
  %   identity and \ indicating the left-multiplication with the inverse.
  % - The left-multiplication (V-I)*x can be calculated efficiently.
  % - The eigenvalues of I-V are similar.
  % Calculations are performed in constant space. During the iteration,
  % space is allocated for two additional copies of the working solution.
  %
  % Multiple problems can be solved in parallel by specifying a b with
  % multiple columns. Run this function without arguments for an example.
  %
  %
  % Input arguments:
  %
  % LpIinv: The matrix inv(L+I) or a function that calculates
  %    the left-multiplication with it in a memory efficient way. The
  %    function should take one input argument and operate on its columns.
  %
  % ImV: The matrix (I-V) or a function that calculates
  %    the left-multiplication with it in a memory efficient way. The
  %    function should take one input argument and operate on its columns.
  %
  % b: The right-hand side of the linear set of equations.
  %
  % tol: The tolerance level. The iteration is terminated when the relative
  %    residue has decreased below this threshold. Default: 1e-6;
  %
  % maxit: The maximum number of iterations before the iteration is
  %    terminated. Each iteration requires to calls to ImV and
  %    one call to LpIinv, as well as 3 additions and the 
  %    calculation of the residual norm.
  %
  % x0: The starting value of the iteration. Default: 0.
  %
  % callback: A function, cont = callback(iter, relres, x, dx), that is 
  %     called after each iteration with the optional input arguments
  %     (iter, relres, x, dx) and a boolean output argument. The iteration
  %     continues as long as the callback function returns 'true'. The
  %     input arguments maxit and tol are ignored when a callback function
  %     is provided.
  %     Example usage:
  %         function cont = check_progress(iter, relres)
  %             if mod(iter, 100) == 0
  %                 fprintf('Iteration %d: error %0.6f%%.\n', [iter, 100*mean(relres(:))]);
  %             end
  %
  %             cont = iter < maxit && any(relres(:) > tol);
  %         end
  %
  %         callback = @check_progress
  %
  % display_progress: A boolean indicating whether the default callback
  %     displays progress every 100 iterations. Default: true
  %
  % Output arguments:
  %
  % x: The result as a column vector, or an nd-array if the right-hand side,
  %    b, is an nd-array.
  %
  % flag: The convergence flag:
  %       0: all converged to the set tolerance,
  %       1: did not reach convergence or was stopped by a user-specified callback,
  %       2: the iteration diverges. Either V is not a contraction, or the
  %          real part of H=L+V is not positive definite.
  %
  % relres: The relative residual, or an nd-array if the right-hand side,
  %    b, is an nd-array.
  %
  % iter: A scalar indicating the number of completed iterations.
  %
  % resvec: A column vector with the absolute residuals after each iteration, or
  %    an nd-array of column vectors when the right-hand side is an nd-array.
  %
  
  asymptotic_stop_criterion_enabled = false;  % Allow it to activate. TODO: fix
  asymptotic_stop_criterion_active = false;  % Start inactive until condition is met
  
  % Default input arguments
  test_vector_length = 1000;
  if nargin < 1 || isempty(LpIinv)
    % Define a test diagonal to convolve with
    ld = (10/test_vector_length) * [1:test_vector_length].';
    ld = (ld - ld(floor(1+end/2))).^2;
    ld = ld ./ ld(1);
%     L_plus_identity = @(x) ifft(fft(x) .* (ld + 1));
    LpIinv = @(y) ifft(fft(y) ./ (ld + 1));
  end
  if nargin < 2 || isempty(ImV)
    vd = (1.0 * rand(test_vector_length, 1) + 0.0) .* exp(1i * pi * (rand(test_vector_length, 1) - 0.5));
    ImV = @(x) x - vd .* x;
  end
  if nargin < 3 || isempty(b)
    b = rand(test_vector_length, 2);
  end
  if nargin < 4 || isempty(tol)
    tol = 2 * eps(class(b));  % 1e-6;
  end
  if nargin < 5 || isempty(maxit)
    maxit = Inf;
  end
  if nargin < 6 || isempty(x0)
    x0 = 0;
  end
  if nargin < 7
    callback = [];
  end
  if nargin < 8 || isempty(display_progress)
    display_progress = true;
  end
  
  % Standarize inputs
  if isnumeric(LpIinv)
    % LpIinv = @(y) (L + eye(size(L))) \ y;
    LpIinv = @(y) LpIinv * y;
  end
  if isnumeric(ImV)
    %ImV = @(x) x - V * x;
    ImV = @(x) ImV * x;
  end
  
  problem_shape = size(b);
  
  function n = norm_single(x)
    n = norm(x);
  end
  function n = norm_generic(x)  % This can be significantly (>10%) slower
    n = zeros(1, problem_shape(2));
    for problem_idx = [1:problem_shape(2)]
      n(problem_idx) = norm(x(:, problem_idx));
    end
  end
  if problem_shape(2) == 1
      column_norm = @norm_single;
  else
      column_norm = @norm_generic;
  end
  
  if any(column_norm(b).^2 < eps(class(b)))
    error('The right hand side (b) contains a vector with a norm close to zero.');
  end
    
  % If possible, avoid initial multiplications with 0
  if all(x0(:) == 0)
    % Short-cut the first iteration
    x = ImV(LpIinv(b));
    residue = column_norm(x);
    iter = 1;
  else
    x = x0;
    clear x0;
    residue = Inf;
    iter = 0;
  end
  
  % Define a default function to check convergence
  function cont = default_callback(iter, relres)
    %
    % The default callback function only takes two arguments
    %
    if display_progress && mod(iter, 100) == 0
      fprintf('Iteration %d: error %0.6f%%.\n', [iter, 100*mean(relres(:))]);
    end
    cont = iter < maxit && any(relres(:) > tol);
  end
  if isempty(callback)
    % Define the callback
    callback = @(iter, relres) default_callback(iter, relres);
  end
  
  % Initialize iteration
  norm_dx = [];
  resvec = zeros([0, problem_shape(2:end)]);
  cont = true;
  while cont  % continue until the callback says otherwise
    % Calculate the correction dx = Gamma(1-V)x - (1-V)x + Gamma b
    dx = ImV(LpIinv(ImV(x) + b) - x);  % Gamma((1-V)x + b) - (1-V)x
    % Update solution  x = Gamma(V-1)x + Vx - Gamma b
    x = x + dx;  % Do not deallocate dx and temporary variable, we will recycle the allocated memory in the next loop
    iter = iter + 1;
    relres = [];
    norm_x = [];
    
    % Convergence check and feedback
    previous_residue = residue;
    
   % Asymptotic stop criterion
   % Insufficiently accurate at early iterations?
    norm_prev_dx = norm_dx;
    norm_dx = column_norm(dx);
    residue = norm_dx;  % Simple stop criterion
    if asymptotic_stop_criterion_enabled
      if ~isempty(norm_prev_dx)
        norm_x = column_norm(x);
        estimate_eig_M = norm_dx ./ norm_prev_dx;
        if ~asymptotic_stop_criterion_active
%           inner_product_update_and_solution = abs(dot(dx, x, 1));
          % angle_between_update_and_solution = acos(inner_product_update_and_solution ./ (norm_dx .* norm_x));
          inner_product_between_updates = abs(dot(dx, prev_dx, 1));
          % angle_between_updates = acos(inner_product_between_updates ./ (norm_dx .* norm_prev_dx));            
%           if all(inner_product_update_and_solution < 0.01 .* norm_dx .* norm_x)
          if all(inner_product_between_updates > 0.9 .* norm_dx .* norm_prev_dx)
            asymptotic_stop_criterion_active = true;
            % logMessage('Estimate of |E_M| = %0.6f, residue = %0.3e, relative residue = %0.3e', [estimate_eig_M, estimate_residue, relres]);
          end
        end
        if asymptotic_stop_criterion_active
          residue = norm_dx ./ (estimate_eig_M.^-1 - 1);
        end
      end
      prev_dx = dx;  % Only needed for the angle calculation above
    end
    
    if nargout >= 5  % If requested, store all the residuals
      resvec(end+1, :) = residue(:, :);
    end
    if all(residue <= previous_residue)
        % Call callback function without redundant arguments
        callback_args = {};
        if nargin(callback) >= 1
            callback_args{end+1} = iter;
        end
        if nargin(callback) >= 2
            if isempty(norm_x)
                norm_x = column_norm(x);
            end
            relres = residue ./ norm_x;
            callback_args{end+1} = relres;
        end
        if nargin(callback) >= 3
            callback_args{end+1} = x;
        end
        if nargin(callback) >= 4
            callback_args{end+1} = dx;
        end
        cont = callback(callback_args{:});
    else
      % Override callback result when divergence detected
      warning('Divergence detected, either V is not a contraction, or the real part of (L+V) is not positive definite!');
      cont = false;
    end
  end  % while
  
  if nargout >= 3 && isempty(relres)
    if isempty(norm_x)
      norm_x = column_norm(x);
    end
    relres = residue ./ norm_x;
  end
  
  % Report convergence as with other Matlab functions
  if any(residue > previous_residue)
    flag = 2;  % divergence detected
  elseif any(relres > tol)
    flag = 1;  % did not converge to tolerance before the maximum number of iterations reached
  else
    flag = 0;  % converged to tolerance
  end
  
end  % function
