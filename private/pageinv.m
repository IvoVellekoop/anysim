function X = pageinv(X)
%PAGEINV(X) Page-wise matrix inversion
%   Each page (i.e. the first two dimensions) of N-dimensional array X
%   is treated as a matrix that is inverted. 
%   SEE ALSO pagemtimes
%
    if isa(X, 'gpuArray') || isa(X, 'distributed')
        X = pagefun(@inv, X);
    else % why is pagefun not implemented for ordinary arrays?
        [U, S, V] = pagesvd(X, 'vector');
        X = pagemtimes(V./pagetranspose(S), 'none', U, 'ctranspose');
    end
end

