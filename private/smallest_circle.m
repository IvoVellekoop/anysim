function [center, radius] = smallest_circle(points, tolerance)
%SMALLEST_CIRCLE Finds the smallest circle that encompasses all points
%   This function uses a very efficient, yet approximate algorithm.
%   The algorithm is guaranteed to return a circle that includes all
%   points. However, this circle may be slightly bigger than optimal,
%   as specified by the optional tolerance parameter (default 1E-10).
%
%   points = points to fit the circle around (real part=x, imaginary part = y) 
%   tolerance = termination tolerance (optional, default 1E-10). The
%               returned radius may be 'tolerance' larger than optimal.
%   returns center point (real part=x, imaginary part = y) and radius of the circle
%
% By default a (slow) reference implementation is used. 
% Build the accelerators (see readme.md) to use a highly optimized version.
%
% Copyright 2019-2022. Ivo Vellekoop, University of Twente
%


% Define end result (smallest circle containing all points) as [c, r]
%
% t_o === [p_o(1:3), r_o] = inflated triangle. 
%  * t_o includes c
%  * r_o >= r
%  * p_o(1:3) are points from the set
%
    % 1. calculate [c_c, r_c], the smallest circle around p_o(1:3), (using a
    %    trivial algorithm)
    %
    % 2. calculate conjugated inflated triangle t_c === [c_c(1:3), r_o]
    %
    % 3. for all points
    %   calculate distances to [center, c_c(1:3)]
    %   if min(|o_c - c(1:3)|) < r_o (i.e. the point is inside t_o)
    %       remove point
    %   else
    %       find points with largest distance for each center 
    %           outliers: [o_c, o(1), o(2), o(3)]
    %   end
    %
    % 3. if |o_c - center| < radius + tolerance, 
    %       return [center, radius] as the result
    %    else
    %       r_o => min |o(1:3) - c(1:3)|
    %       p_o => o(1:3)
    %    end
    % (repeat from 1) 
    points = points(:);
    
    % special case for points on the real line to ensure returned center is
    % real (it could otherwise have a small imaginary component due to rounding errors)
    %    
    if isreal(points)
        pmin = min(points);
        pmax = max(points);
        center = (pmin + pmax) / 2;
        radius = pmax - center;
        return;
    end


    N_reads = 0; % number of elements read from memory
    if nargin < 2
        tolerance = 1E-10;
    end
    draw = false;

% Step 0, pick four initial corner points based on bounding box
    corner_i = zeros(4,1);
    [~, corner_i(1)] = min(real(points));
    [~, corner_i(2)] = max(real(points));
    [~, corner_i(3)] = min(imag(points));
    [~, corner_i(4)] = max(imag(points));
    p_o = points(corner_i);
    width = real(p_o(2) - p_o(1));
    height = imag(p_o(4) - p_o(3));
    r_o = sqrt(width^2 + height^2) / 2;
    center = 0.5 * (real(p_o(1) + p_o(2)) + 1i * imag(p_o(3) + p_o(4)));
    N_reads = N_reads + numel(points);
    
    for it=1:50
        % step 1
        % here, p_o contains up to 7 points, pick the 2-3 that correspond
        % to the smallest circle enclosing them all
        % sort in order of increasing distance from center since it is
        % more likely that the final circle will be built from the points
        % further away.
        %
        [~, ind] = sort(abs(p_o-center));
        [center, radius, p_o] = smallest_circle_brute_force(p_o(ind));
        
        % step 2
        c_c = conjugate_inflated_triangle(p_o, r_o);
        
        % 2a: select points
        distances = abs(points - [center; c_c].');
        keep = max(distances, [], 2) > r_o;% - tolerance;
        N_reads = N_reads + numel(keep);
        
        % 2b: determine outliers
        [r_out, outliers_i] = max(distances, [], 1);
        if r_out(1) < radius + tolerance
            radius = r_out(1);
            return;
        end
        
        %%
        if draw
            plot(points, '.r'); hold on;
            plot(points(keep), '.g');
            plot(center, 'o');
            plot(p_o, '+');
            plot(c_c, 'd'); 
            draw_circle(p_o(1), r_o);
            draw_circle(p_o(2), r_o);
            if length(p_o) > 2
                draw_circle(p_o(3), r_o);
            end
            plot(points(outliers_i), '*');
            plot(5+3i, 's'); %exact center
            axis image;
            hold off;
        end
        outliers = points(outliers_i);
        r_o = min(min(r_out), r_o);
        points = [points(keep); p_o];
        p_o = [outliers; p_o];        
    end
    warning('Did not converge in 50 iterations, tolerance is set too low for machine precision?');
end

function c_c = conjugate_inflated_triangle(points, r)
    c_c = zeros(3,1);

    Np = length(points);
    if Np == 2
        B = points(1);
        C = points(2);
        M_start = 0.5 * (C+B);
        M_dir = 1i * (C-B) / abs(C-B);
        w = abs(C - M_start);
        alpha = sqrt(r^2 - w^2);
        c_c(1) = M_start - alpha * M_dir;
        c_c(2) = M_start + alpha * M_dir;
        
        %not needed, but we can pick one point 'for free'
        M_dir = (C-B) / abs(C-B);
        c_c(3) = B + M_dir * r;
        return;
    end
    
    ss = signed_surface(points);
    if ss < 0 % reverse direction of circle (otherwise we get the solutions outside of the circle)
        tmp = points(1);
        points(1) = points(3);
        points(3) = tmp;
    end
    
    for p=1:Np
        % For conjugating point A, define mid-point M_CB and line M_CB-A
        % c_A is on this line, at a distance of r from point B (and point
        % C)
        %
        % c_A = M_CB + alpha * (A - M_CB) / |A-M_CB|  with alpha >= 0 
        % w = |C - M_CB|
        % h = alpha
        % w^2 + alpha^2 = r^2
        % alpha = sqrt(r^2 - w^2)
        %
        B = points(p);
        C = points(mod(p, Np) + 1);
        M_start = 0.5 * (C+B);
        M_dir = 1i * (C-B) / abs(C-B);
        w = abs(C - M_start);
        alpha = sqrt(r^2 - w^2);
        c_c(p) = M_start + alpha * M_dir;
    end
end

function S = signed_surface(points)
    S = imag((points(1) - points(2)) .* conj(points(3)-points(2))); 
end

function draw_circle(center, radius, varargin)
    rectangle('Position', [real(center)-radius, imag(center)-radius, 2*radius, 2*radius], 'Curvature', 1, varargin{:});
end