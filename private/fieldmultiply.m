function x = fieldmultiply(A, x)
%FIELDMULTIPLY(A, x) Computes A x
%   There are three possible cases:
%       1. A and x are both scalar fields.
%       2. A is a field of diagonal matrices, and x is a vector field.
%       3. A is a field of full matrices, and x is a vector field.
%   
%   To detect which case is applicable:
%   * if A and x have the same dimensions we have case 1 or 2
%   * if A is a scalar or a vector, we have case 1 or 2
%   * in any other situation we have case 3
%
% Note that there is one edge case: suppose we have a 3-element vector
% field of size 3 x 1, and A is a full matrix. Then both x and A have size 3x3
% and we incorrectly treat this as case 1 or 2.
% We cannot distinguish the situations because MATLAB removes trailing
% singleton dimensions.
% To prevent this situation, a check is done in the constructor of the
% simulation object to prevent the latter case (i.e. x can never be a
% square matrix if it holds a vector field)

    if isequal(size(A), size(x)) || size(A, 2) == 1
        x = A .* x;
    else
        sz = size(x);
        x = reshape(x, [sz(1) 1 sz(2:end)]);
        x = pagemtimes(A, x);
        x = reshape(x, sz);
    end
end

