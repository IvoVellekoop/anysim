function N = pagenorm(X)
%PAGEMNORM(X) Computes the maximum induced 2-norm for all pages in X
%   Each page (i.e. the first two dimensions) of N-dimensional array X
%   is treated as a matrix for which the norm is computed. The largest
%   norm that is found is returned.
%   SEE ALSO pagemtimes, pageinv
%
    if isa(X, 'distributed')
        N = max(pagefun(@norm, X), [], 'all');
    else
        N = max(pagesvd(gather(X), 'vector'), [], 'all');
    end
end

