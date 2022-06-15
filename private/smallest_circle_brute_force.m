function [center, radius, corners] = smallest_circle_brute_force(points, tolerance)
% Determines the smallest circle encompassing all points
%   Points are specified as x+1i*y
%   Returns center point, radius, and values of the 1-3 edge points
%   that define the final circle
%
% Only used in reference implementation and test. Not used when
% optimized accelerator is used.
% Note this algorithm is very inefficient for more than 4 points
    if nargin < 2
        tolerance = 1E-10;
    end

% special cases for one and two points
    N = length(points);
    if N == 1
        center = points;
        radius = 0;
        corners = points;
        return;
    elseif N == 2
        center = 0.5 * (points(1) + points(2));
        radius = abs(points(1)-center);
        corners = points;
        return;        
    end

% Remove one point and recursively construct a smallest circle for the
% remaining points. If the removed point is inside that circle, return
% if the removed point is not in the circle, repeat with a different point
% omitted. First check if it is possible to construct a circle from just 2 points,
% including the third
% todo: faster check to see if two or three points are needed?
    Ns = 1:N;
    for p=1:N
        reduced = points(Ns~=p); %profiler shows this is actually the slowest line in this function?!
        [center, radius, corners] = smallest_circle_brute_force(reduced, tolerance);
        if abs(points(p)-center) <= radius + tolerance
            return; % The omitted point is inside the circle spanned by the other points
        end
    end
    
    % if we get here, no suitable subset was found. This is only possible
    % for 3 points

    %% All three points are edge points
%(squared) distance from center to all three points is equal:
%    cA2 = |c-A|^2 = |c|^2 + 2 Re (A^* c) + |A|^2
%    cB2 = |c-B|^2 = |c|^2 + 2 Re (B^* c) + |B|^2
%    cC2 = |c-C|^2 = |c|^2 + 2 Re (C^* c) + |C|^2
    
% from cA2 - cB2
%    2 Re (A^* c) - 2 Re (B^* c) = |A|^2 - |B|^2
% from cA2 - cC2
%    2 Re (C^* c) - 2 Re (C^* c) = |A|^2 - |C|^2
%
% write real/imaginary parts
%
%   2 (Re(A)-Re(B)) * Re(c) + 2 (Im(A)-Im(B) * Im(c)) = |A|^2 - |B|^2
%   2 (Re(A)-Re(C)) * Re(c) + 2 (Im(A)-Im(C) * Im(c)) = |A|^2 - |C|^2
%
% now write in matrix form and let MATLAB solve it
    A = points(1);
    B = points(2);
    C = points(3);
    M = 2*[real(A)-real(B), imag(A)-imag(B); real(A)-real(C), imag(A)-imag(C)];
    b = [abs(A)^2 - abs(B)^2; abs(A)^2 - abs(C)^2];
    c = M\b;
    center = c(1) + 1i * c(2);
    radius = abs(A-center);
    corners = points;
end