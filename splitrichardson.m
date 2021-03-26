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
  % resvec: A column vector with relative residual after each iteration, or
  %    an nd-array of column vectors when the right-hand side is an
  %    nd-array.
  %
  
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
  if any(norm2(b) < eps(class(b)))
    error('The right hand side (b) contains a vector with a norm close to zero.');
  end
    
  % Define a function that calculates the relative norm column-by-column
  function n = rel_norm_single(dx, x)
    n = norm(dx) ./ norm(x);
  end
  function n = rel_norm_generic(dx, x)  % This can be significantly (>10%) slower
    n = zeros(1, problem_shape(2));
    for problem_idx = [1:problem_shape(2)]
      n(problem_idx) = norm(dx(:, problem_idx)) ./ norm(x(:, problem_idx));
    end
  end
  if problem_shape(2) == 1
      rel_norm = @rel_norm_single;
  else
      rel_norm = @rel_norm_generic;
  end

  % If possible, avoid initial multiplications with 0
  if all(x0(:) == 0)
    % Short-cut the first iteration
    x = ImV(LpIinv(b));
    relres = 1;
    iter = 1;
  else
    x = x0;
    clear x0;
    relres = Inf;
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
  dx = [];
  norm_dx = [];
  resvec = zeros([0, problem_shape(2:end)]);
  cont = true;
  while cont  % continue until the callback says otherwise
%     prev_dx = dx;  % Only needed for the asymptotic stop criterion below
    % Calculate the correction dx = Gamma(1-V)x - (1-V)x + Gamma b
    dx = ImV(LpIinv(ImV(x) + b) - x);  % Gamma((1-V)x + b) - (1-V)x
    % Update solution  x = Gamma(V-1)x + Vx - Gamma b
    x = x + dx;  % Do not deallocate dx and temporary variable, we will recycle the allocated memory in the next loop
    iter = iter + 1;
    
    % Convergence check and feedback
    previous_relres = relres;
    relres = rel_norm(dx, x);
    
%    % Asymptotic stop criterion
%    % Insufficiently accurate at early iterations?
%     norm_prev_dx = norm_dx;
%     norm_dx = sqrt(norm2(dx));
%     if ~isempty(norm_prev_dx)
%         estimate_eig_M = norm_dx ./ norm_prev_dx;
%         % norm_x = sqrt(norm2(x));
%         % angle_between_update_and_solution = acos(abs(sum(conj(dx) .* x, 1)) ./ (norm_dx .* norm_x));
%         angle_between_updates = acos(abs(sum(conj(dx) .* prev_dx, 1)) ./ (norm_dx .* norm_prev_dx));
%         % logMessage('%0.1f degrees between update and solution. %0.1f degrees between current update and previous update.', [angle_between_update_and_solution, angle_between_updates] .* 180 ./ pi);
%         estimate_residue = estimate_eig_M .* (1 - estimate_eig_M) .* norm_dx;
%         norm_x = sqrt(norm2(x));
%         relres_actual = estimate_residue ./ norm_x;
%         % logMessage('Estimate of |E_M| = %0.6f, residue = %0.3e, relative residue = %0.3e', [estimate_eig_M, estimate_residue, relres]);
%         if all(angle_between_updates < 0.1)
%             % logMessage('Using asymptotic stop criterion.');
%             relres = relres_actual;
%         end
%     end
    
    if nargout >= 5 
      resvec(end+1, :) = relres(:, :);
    end
    % Call callback function without redundant arguments
    switch nargin(callback)
      case 0
        cont = callback();
      case 1
        cont = callback(iter);
      case 2
        cont = callback(iter, relres);
      case 3
        cont = callback(iter, relres, x);
      otherwise
        cont = callback(iter, relres, x, dx);
    end  % switch
    % Override callback when divergence detected
    if any(relres > previous_relres)
      warning('Divergence detected, either V is not a contraction, or the real part of (L+V) is not positive definite!');
      cont = false;
    end
  end  % while
  
  % Report convergence as with other Matlab functions
  if any(relres > previous_relres)
    flag = 2;  % divergence detected
  elseif any(relres > tol)
    flag = 1;  % did not converge to tolerance before the maximum number of iterations reached
  else
    flag = 0;  % converged to tolerance
  end
  
  % The resvec return argument contains absolute residues
  if nargout >= 5 
    resvec = resvec .* sqrt(norm2(x));
  end
  
end  % function

function n2 = norm2(x)
  %
  % Calculate the square of the l2-vector norm in parallel for all columns
  %
  n2 = sum(abs(x).^2, 1);
end
