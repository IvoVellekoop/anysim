function [center, radius] = smallest_circle(points)
%SMALLEST_CIRCLE Finds the smallest circle that encompasses all points
%   At the moment, this function only uses a very crude approximation:
%   just take center as the point exactly between the minimum and maximum 
%   value in both the x(real) and y(imaginary) direction.
%
    x_min = min(real(points));
    y_min = min(imag(points));
    x_max = max(real(points));
    y_max = max(imag(points));
    center = 0.5 * (x_min+x_max) + 0.5i * (y_min+y_max);
    radius = max(abs(points-center));
end