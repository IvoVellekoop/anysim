function mustBeSquare(A)
%MUSTBESQUARE Validate that A is a square matrix (or a scalar)
%   MUSTBESQUARE(A) throws an error if A is not square.
%
%   Class support:
%   All numeric classes, logical
%   MATLAB classes that define these methods:
%       gt, isreal, isnumeric, islogical
%
    if size(A, 1) ~= size(A, 2) || ~ismatrix(A)
        throwAsCaller(MException('ANYSIM:validators:mustBeNumericSquare', "Matrix must be square"));
    end
end
