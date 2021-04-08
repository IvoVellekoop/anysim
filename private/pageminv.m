function X = pageminv(X)
%PAGEMINV(X) Page-wise matrix inversion
%   Each page (i.e. the first two dimensions) of N-dimensional array X
%   is treated as a matrix that is inverted. 
%   SEE ALSO pagemtimes
%
    if isa(X, 'gpuArray') || isa(X, 'distributed')
        X = pagefun(@inv, X);
    else % why is pagefun not implemented for ordinary arrays?
        sX = size(X);
        X = X(:,:,:);
        for n=1:size(X, 3)
            X(:,:,n) = inv(X(:,:,n));
        end
        X = reshape(X, sX);
    end
end

