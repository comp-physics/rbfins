function d = drounded_rectangle(p, x_center, y_center, width, height, radius)
%DROUNDED_RECTANGLE Distance function for rounded rectangle optimized for DistMesh
%
% This function computes the signed distance from points to a rounded
% rectangle obstacle using a smooth, well-conditioned formulation that
% works reliably with DistMesh.
%
% INPUTS:
%   p        - Points [N x 2] where distance is evaluated
%   x_center - X-coordinate of rectangle center
%   y_center - Y-coordinate of rectangle center
%   width    - Rectangle width (x-direction)
%   height   - Rectangle height (y-direction)
%   radius   - Corner rounding radius
%
% OUTPUTS:
%   d        - Signed distance values [N x 1]

% Translate to origin-centered rectangle
px = p(:, 1) - x_center;
py = p(:, 2) - y_center;

% Half dimensions
a = width / 2;
b = height / 2;

% Clamp radius to prevent invalid geometry
max_radius = min(a, b) * 0.9; % Leave some margin
radius = min(radius, max_radius);

% Compute distance using smooth formulation
% Distance to rectangle edges (before rounding)
dx = max(abs(px) - (a - radius), 0);
dy = max(abs(py) - (b - radius), 0);

% Smooth distance function for rounded rectangle
d = sqrt(dx.^2 + dy.^2) - radius;

% Handle interior points correctly
interior_mask = (abs(px) <= a - radius) | (abs(py) <= b - radius);
d(interior_mask) = -radius + max(abs(px(interior_mask)) - (a - radius), ...
                                 abs(py(interior_mask)) - (b - radius));

% Final adjustment for interior
inside_box = (abs(px) <= a) & (abs(py) <= b);
d(inside_box) = min(d(inside_box), -min(a - abs(px(inside_box)), b - abs(py(inside_box))));

end
