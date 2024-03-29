classdef (Abstract) AnySim
    %ANYSIM Simulation system for linear operator equations
    %   (c) 2019. Ivo Vellekoop
    %
    % This is the base class for all simulation objects. It implements
    % the split-Richardson algorithm in an abstract way.
    % It is the responsibility of the constructor of the simulation
    % object to implement the operators for the medium,
    % propagator, source, transform and scaling
    %
    properties (SetAccess = protected)
        % operators (see readme.md) for a detailed explanation
        medium      % operator V+1. This object implements two functions:
                    % mix(phi, prop_phi) = (V+1)^2 prop_phi + (V+1) phi
                    % mix_final(phi, prop_phi) = (V+1) prop_phi + phi
                    % also includes fields for preconditioning matrices Tl, V0 and Tr
        propagator  % operator (L-1)^-1
        medium_adj = @AnySim.medium_adj_err
        propagator_adj = @AnySim.propagator_adj_err
        L           % Function handle or matrix for the 
        L_adj       % 'forward' operator L.
                    % Since this operators are not needed by anysim
                    % itself, it is only generated when 
                    % opt.forward_operator == true
        Tl; Tr      % pre-preconditioning operators Tl, Tr
        V0          % centering potential
        opt         % options
    end
        
    methods
        function obj = AnySim(opt)
            arguments
                opt (1,1) AnySimOptions
            end
            obj.opt = opt;
        end
        
        function x = preconditioner(obj, x, state)
            % SIM.PRECONDITIONER(X) applies the selected (inverse) preconditioner on X
            % By default, this is B(L+1)^(-1) X
            %
            % Note: this function is not used by the anysim
            % algorithm itself, but it can be used to compare
            % anysim so other algoriths (such as GMRES)
            %
            switch obj.opt.preconditioner
            case "none"
                x = x;
                if (norm(x(:)) < 1E-4)
                    warning("near-zero");
                end
            case "shift" 
                % Bai, Zhong-zhi, et al. “A SHIFT-SPLITTING PRECONDITIONER FOR NON-HERMITIAN POSITIVE DEFINITE MATRICES.” Journal of Computational Mathematics, vol. 24, no. 4, 2006, pp. 539–52.
                % Γ = 1/2 (A + αI)
                % Γ^-1 A = 2 (A + αI)^-1 (A + αI - αI) = 2 - 2α(A + αI)^-1  
                % The preconditioner is evaluated with a nested iterative algorithm.
                % However, a preconditioner is needed for this nested algorithm too.
                % For that, we use the AnySim preconditioner (since that is
                % the only one we know that guarantees convergence under
                % all circumstances)
                % therefore, here we solve:
                % B (L+1)^-1 y = B (L+1)^-1 1/2 (A + αI) x
                %xcopy = x;
                x = obj.medium(obj.propagator(x));
                [x,flag] = bicgstab(@(u) shift(obj, u, state), x(:), 1E-4, 1000);
                if flag ~= 0
                    if flag == 4 && ~state.internal_iteration_failed
                        state.internal_iteration_failed = true;
                        warning(["Inner iteration converged to " norm(x(:))]);
                    else
                        error("Inner iteration did not converge");
                    end
                end
                x = 2 * reshape(x, [obj.grid.N_u, 1]);
                %xback = (obj.L(x) + (1+obj.opt.preconditioner_shift) * x - obj.medium(x))/2;
            case "hermitian"
                error ("not implemented")
            case "skew-hermitian"
                error ("not implemented")
            case "moborn"
                x = obj.propagator(x);
                x = obj.medium(x);
            case "circulant"
                x = obj.propagator(x);
            end
        end
        function t1 = shift(obj, u, state)
                % compute preconditioned shifted forward operator:
                % B (L+1)^-1(A + α) = B (L+1)^-1 (L+1+α-B) =
                %   B[1 + (L+1)^-1 (α-B)]
                %
                % 
                u = reshape(u, [obj.grid.N_u, 1]);
                t1 = obj.medium(u + obj.propagator(obj.opt.preconditioner_shift * u - obj.medium(u)));
                %t1 = obj.medium(obj.propagator(obj.forward(u) + obj.opt.preconditioner_shift * u));
                t1 = t1(:);
                state.internal_iteration = state.internal_iteration + 1;
               % state.next(u(:), t1(:)); % keep track of evaluations of operator L
               % state.next(u(:), t1(:)); % keep track of evaluations of operator L
        end        
        
        function [u, state] = exec(obj, b)
            arguments
                obj AnySim
                b
            end
            % EXEC executes the simulation for the given source
            % See README.md for a detailed explanation of the algorithm

            %% Initialize state and source
            % The iteration is implemented as:
            % t1 = G u + s            Medium.mix_source
            % t1 -> Li t1             Propagator.propagate
            % u -> u + G (t1 - u)     Medium.mix_field
            [u, state] = obj.start();
            
            while state.running
                % t1 => B u + b
                t1 = obj.medium(u) + b; 
                
                % t1 => (L+1)^-1 t1
                t1 = obj.propagator(t1);
                
                % u + B (t1-u)
                t1 = obj.medium(u - t1); % residual
                state.next(u, t1);
                u = u - obj.opt.alpha * t1;
            end
            
            % u -> Tr u (convert to non-scaled solution)
            % + optional final processing (in grid-based simulations the solution
            % is cropped to the roi)
            u = obj.finalize(u);
            state.finalize();
        end

        function t1 = forward(obj, u, state)
            %% SIM.FORWARD(U) applies the forward operator A to U
            t1 = obj.L(u) + u - obj.medium(u); % L + V = L + 1 - (1-V)
            if nargin > 2
                state.next(u, t1);
            end
        end

        function [f, state] = preconditioned(obj)
            % SIM.PRECONDITIONED(U) returns a function to evaluate the preconditioned operators
            % 
            % For the AnySim proeconditioner, evaluation of the preconditioned operator is optimized
            % so that the forward operator is eliminated. For other
            % preconditioners, this function first applies the forward
            % operator and then the preconditioner separately (see
            % obj.preconditioner)

            % Usage:
            % A = sim.preconditioned
            % A(u)                      % computes (1-V)(L+1)^(-1) (L+V) u
            % 
            [~, state] = obj.start();
            switch obj.opt.preconditioner
            case "moborn"
                f = @(u) moborn(obj, state, u);
            case "shift"
                f = @(u) shift_preconditioned(obj, state, u);
            otherwise
                f = @(u) preconditioned(obj, state, u);
                %case "circulant"
                % Γ^-1 A = (L+1)^-1 A  
                %case "hermitian"
                %case "skew-hermitian"
            end


            function t1 = preconditioned(obj, state, u)
                t1 = obj.preconditioner(forward(obj, u, state), state);
            end

            function t1 = moborn(obj, state, u)
                t1 = obj.medium(u);         % (1-V)u
                t1 = obj.propagator(t1);    % (L+1)^(-1) (1-V)u
                t1 = obj.medium(u - t1);    % (1-V) (u-t1)
                state.next(u, t1);
            end
            function t1 = shift_preconditioned(obj, state, u)
                % Γ^-1 A = 2 (A + αI)^-1 A = 2 (A + αI)^-1 (A + αI - αI) 
                % Γ^-1 A = 2 - 2α(A + αI)^-1 
                % to solve the inverse, we need to precondition the
                % A + αI operator
                x = obj.medium(obj.propagator(u));
                [x,flag] = bicgstab(@(u) shift(obj, u, state), x(:), 1E-4, 1000);
                if flag ~= 0
                    if flag == 4 && ~state.internal_iteration_failed
                        state.internal_iteration_failed = true;
                        warning("Inner iteration did not converge");
                    else
                        error("Inner iteration did not converge");
                    end
                end
                t1 = 2 * u - 2 * obj.opt.preconditioner_shift * reshape(x, [obj.grid.N_u, 1]);
                state.next(u, t1);
            end
        end
    end
    methods (Abstract)
        u = finalize(obj, u)
    end
    methods (Abstract, Access=protected)
        [u, state] = start(obj)
    end
    methods (Static, Access=private)
        function u = medium_adj_err(u) %#ok<INUSD> 
            error("B^* (adjoint medium operator) is not defined");
        end
        function u = propagator_adj_err(u) %#ok<INUSD> 
            error("(L^*+1)^-1 (adjoint propagetor) is not defined");
        end
    end
end


