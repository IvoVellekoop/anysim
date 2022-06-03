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
                    % 'forward' operator L.
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
        
        function x = preconditioner(obj, x)
            % SIM.PRECONDITIONER(X) returns B(L+1)^(-1) X
            %
            % Note: this function is not used by the anysim
            % algorithm itself, but it can be used to compare
            % anysim so other algoriths (such as GMRES)
            %
            x = obj.propagator(x);
            x = obj.medium(x);
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
        
        function [f, state] = preconditioned(obj)
            % SIM.PRECONDITIONED(U) returns a function to evaluate the preconditioned operators
            % 
            % Usage:
            % A = sim.preconditioned
            % A(u)                      % computes (1-V)(L+1)^(-1) (L+V) u
            % 
            % Note: (L+1)^(-1) L = 1-(L+1)^(-1)
            % So: (1-V)(L+1)^(-1) (L+V)  
            %  =  (1-V)[1 - (L+1)^(-1)] + (1-V)(L+1)^(-1) V 
            %  =  (1-V)[1-(L+1)^(-1)(1-V)]
            %  = B[1-(L+1)^{-1} B]
            % 
            % and the adjoint:
            %  =  [1-(1-V*)(L*+1)^(-1)](1-V*)
            %  =  (1-V*)[1-(L*+1)^(-1)](1-V*)
            % (so we can just take the adjoint of V and L and perform the
            % same sequence of operations)
            %
            [~, state] = obj.start();
            f = @(varargin) apply_preconditioned(obj, state, varargin{:}); 

            function t1 = apply_preconditioned(obj, state, u, transpose_flag)
                if nargin == 4 && strcmp(transpose_flag,'transp')
                    t1 = obj.medium_adj(u);         % (1-V*)u
                    t1 = obj.propagator_adj(t1);    % (L*+1)^(-1) (1-V*)u
                    t1 = obj.medium_adj(u - t1);    % (1-V*) (u-t1)
                else
                    t1 = obj.medium(u);         % (1-V)u
                    t1 = obj.propagator(t1);    % (L+1)^(-1) (1-V)u
                    t1 = obj.medium(u - t1);    % (1-V) (u-t1)
                end
                state.next(u, t1);
            end
        end
        
        function [f, state] = operator(obj)
            % SIM.OPERATOR(U) Returns (L+V)U
            %
            % U should be in the domain of V (typically real space)
            %
            % Note: this function is not used by the anysim
            % algorithm itself, but it can be used to compare
            % anysim so other algoriths (such as GMRES), or to
            % compute the final residual ‖(L+V)U - S‖
            %
            % Note: because this operator is usually not needed and 
            % construction may be expensive, it is only constructed when
            % the option .forward_operator == true
            %
            % Note: These are the scaled L and V, without
            % preconditioner
            %            
            if isempty(obj.L) 
                error('No forward operator was generated, set opt.forward_operator=true and verify that this simulation supports forward operator generation');
            end

            [~, state] = obj.start();
            f = @(varargin) forward_operator(obj, state, varargin{:});

            function t1 = forward_operator(obj, state, u, transpose_flag)
                t1 = obj.L(u) + u - obj.medium(u); % L + V = L + 1 - (1-V)
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


