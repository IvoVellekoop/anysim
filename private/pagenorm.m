function X = pagenorm(X)
%PAGEMNORM(X) Page-wise matrix norm
%   Each page (i.e. the first two dimensions) of N-dimensional array X
%   is treated as a matrix that is inverted. 
%   SEE ALSO pagemtimes
%
    if isa(X, 'distributed')
        X = pagefun(@norm, X);
    else % why is pagefun not implemented for ordinary arrays?
        sX = size(X);
        X = X(:,:,:);
        for n=1:size(X, 3)
            X(:,:,n) = norm(X(:,:,n));
        end
        X = reshape(X, sX);
    end
end

