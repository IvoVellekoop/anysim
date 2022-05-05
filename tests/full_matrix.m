function M = full_matrix(A, N)
%FULL_MATRIX(A, N)
%   Converts operator A to an N x N matrix. A should be a function taking
%   a single column vector as input.
    shape = N;
    Nf = prod(N);
    M = zeros(Nf, Nf);
    b = zeros(Nf, 1);
    b(1) = 1;
    for n=1:Nf
        M(:,n) = reshape(A(reshape(b, shape)), [], 1);
        b = circshift(b, [1,0]);
    end
end

