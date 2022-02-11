function [center, radius] = smallest_circle(points)
%SMALLEST_CIRCLE Finds the smallest circle that encompasses all points
%   At the moment, this function only uses a very crude approximation:
%   just take center as the point exactly between the minimum and maximum 
%   value in both the x(real) and y(imaginary) direction.
% (c) Tom Vettenburg
    p = double(points(:));
    center0 = mean(p);
    function r = calc_radius(c)
        r = max(abs(p - c));
    end
    center = fminsearch(@calc_radius, center0);
    radius = calc_radius(center);
end
