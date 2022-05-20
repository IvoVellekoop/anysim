function [center, radius] = smallest_circle(points)
%SMALLEST_CIRCLE Finds the smallest circle that encompasses all points
% (c) Tom Vettenburg
    p = double(points(:));
    center0 = mean(p);
    function r = calc_radius(c)
        cc = c(1) + 1i * c(2);
        r = max(abs(p - cc));
    end
    center = fminsearch(@calc_radius, [real(center0), imag(center0)]);
    radius = calc_radius(center);
    center = center(1) + 1i * center(2);
end
