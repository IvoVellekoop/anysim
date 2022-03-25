classdef MatrixSolve < AnySim
    %MATRIXSOLVE General matrix solver for large matrices
    %   Demonstrates that AnySim also works for non-accretive systems
    %   The system is first made positive-definite by computing A*A'
    %   and then solved with the usual preconditioner.
    %   Note: this code is optimized for sparse matrices. The method
    %   is particularly efficient if the condition number of the matrix is
    %   low.
    %   (c) 2022. Ivo Vellekoop
    properties 
        A
    end
    methods
        function obj = MatrixSolve(A, opt)
            arguments
                % Dimensions of the grid in voxels. When empty, this should be 
                % determined automatically from the passed potential array.
                A (:,:)
                opt (1,1) AnySimOptions = AnySimOptions()
            end
            obj@AnySim(opt);
            obj.A = A;

            % construct operator A * A' - diag(A * A')
            % A is a sparse matrix, we don't want to compute A*A' directly
            % because it reduces the sparsity of the result a lot
            L = sum(abs(A.^2), 2); % = diag (A * A')
            
            %todo: write norm_est function that takes function handle
            %V = @(u) A * (A' * u) - L .* u;
            V = A * A' - diag(L);
            norm_est = normest(V, 1e-2); 
            
            scale = opt.V_max / norm_est;
            obj.Tl = 1;
            obj.Tr = scale; 

            % optimize for sparse matrices:
            A = A * sqrt(scale);
            Adiag = 1 + L * scale;
            
%            B = speye(size(A, 1)) - V * scale;
            Linv = 1./(L * scale + 1);
            
%            obj.medium = @(u) B * u;
            obj.medium = @(u) Adiag .* u - A * (A' * u);
            obj.propagator = @(u) Linv .* u;
        end
        function u = finalize(obj, u)
            u = obj.A' * obj.Tr * u;
        end
    end
    methods (Access = protected)
        function [u, state] = start(obj)
            u = zeros(size(obj.A, 1), 1);
            state = State(obj, obj.opt);
        end
    end
end
