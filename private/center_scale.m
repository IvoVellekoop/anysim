function [Tl, Tr, V0, V] = center_scale(Vraw, Vrmin, Vmax)
%CENTER_SCALE Computes scaling matrices TL and TR, and offset V0 so that
%   ||V|| <= 1, with || || the operator norm and V = TL (VRAW - V0) TR
%   under the constraint that the real(V) >= VMIN for each element in Vmin.
%   V0 is chosen so that the scaling matrices can be 'as large as possible'
%   
%
%   [TL, TR, V0] CENTER_SCALE(VRAW, VMIN) with VMIN a scalar treats VRAW as a scalar field
%           and computes TR and V) such that TR * max(abs(VRAW-V0), [], 'all') = 1.
%           where TR is a positive real scalar and V0 is a complex scalar.
%           TL is always 1.
%
%   [TL, TR, V0] CENTER_SCALE(VRAW, VMIN) with VMIN a column vector treats VRAW as a field of
%           diagonal tensors. The diagonal elements are stored as columns
%           in VRAW. Computes positive real vector TR and complex vector V0 such
%           that max(abs((V-V0) * diag(TR))), [], 'all') = 1, where * denotes
%           vector-matrix multiplication for each individual column of V.
%           TL is always 1.
%
%   [TL, TR, V0] CENTER_SCALE(VRAW, VMIN) with VMIN a square matrix treats VRAW as a field of
%           2-dimensional tensors, with the tensors stored as the first
%           two dimensions ('pages') in V.
%           Computes matrices TL, TR and V0 such that 
%           max(norm(TL * (V-V0) * TR)), [], 'all') = 1, where * denotes
%           matrix-matrix multiplication for each matrix of V, and norm computes
%           the L2 operator norm for each matrix.
%           Note: this uses a heuristic that does guarantee optimality.
%           Note: at the moment, TL and TR are both diagonal.
%

% reshapes V so that the first 2 dimensions always hold a
% scalar/vector/matrix and dimensions 2,3,... are the spatiotemporal
% ones
    dim = sum(size(Vrmin) > 1);
            
    Vreshape = shiftdim(Vraw, dim-2); 
    Vrmin = shiftdim(Vrmin, dim-2); 
    N = size(Vreshape, 1);
    M = size(Vreshape, 2);
    centers = zeros(N, M);
    radii = zeros(N, M);
    for n=1:N
        for m=1:M
            [c, r] = smallest_circle(Vreshape(n, m, :));
            % adjust centers and radii so that the real part of centers+radii >=
            % Vmin
            re_diff = real(c) + r - Vrmin(n, m); 
            if re_diff < 0
                c = c - re_diff/2;
                r = r - re_diff/2;
                warning('slowing down simulation to accomodate boundary conditions');
            end
            centers(n,m) = c;
            radii(n,m) = r;
        end
    end
    
    switch dim
        case 0
            Tl = 1;
            Tr = Vmax/radii;
            V0 = centers;
            V = Tr * (Vraw - V0);
        case 1
            if any(radii < abs(centers) * 1E-6)
                radii = max(radii, abs(c) * 1E-6);
                warning('One of the components of the potential is (near-)constant, using threshold to avoid divergence in Tr');
            end
            Tl = eye(M);
            Tr = diag(Vmax./radii(:));
            V0 = centers(:);
            V = diag(Tr) .* (Vraw - V0);
        case 2
            % check if matrix is near-singular
            % and add a small offset to the singular values if it is
            % warning('untested code!');
            [U,S,V] = svd(radii);
            cS = diag(diag(U' * centers * V));
            if any(S) < abs(cS) * 1E-6
                S = max(S, abs(cS) * 1E-6);
                radii = U*S*V';
                warning('One of the components of the potential is (near-)constant, using threshold to avoid divergence in Tr');
            end
            [P,R,C] = equilibrate(radii);
            Tl = P'*R*P;
            Tr = C;
            V0 = centers;
            V = pagemtimes(pagemtimes(Tl, Vraw - V0), Tr);
                
            % The procedure above does not guarantee that ||V|| = Vmax
            % In this final step, we compute the norm of V everywhere
            % and adjust the global scaling as needed
            Tr = Tr / max(pagenorm(V), [], 'all') * Vmax;
            V = pagemtimes(pagemtimes(Tl, Vraw - V0), Tr);
    end
end

