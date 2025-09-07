function d = dellipse(p, xc, yc, a, b)
%DELLIPSE Distance function for ellipse
%
% This function computes the signed distance from points to an ellipse.
% Negative distances indicate points inside the ellipse.
% Uses a more accurate distance approximation for better numerical stability.
%
% INPUTS:
%   p  - Points [x, y] as an Nx2 matrix
%   xc - x-coordinate of ellipse center
%   yc - y-coordinate of ellipse center
%   a  - Semi-major axis (x-direction)
%   b  - Semi-minor axis (y-direction)
%
% OUTPUTS:
%   d  - Signed distance to ellipse boundary (negative inside)

% Translate points to ellipse-centered coordinates
x = p(:, 1) - xc;
y = p(:, 2) - yc;

% Normalize coordinates by ellipse axes
u = x / a;
v = y / b;

% Compute ellipse function value: u^2 + v^2 for the normalized ellipse
ellipse_val = sqrt(u.^2+v.^2);

% For better numerical stability, use a more accurate distance approximation
% that accounts for the ellipse curvature

% Compute the characteristic length scale at each point
% This varies around the ellipse perimeter
theta = atan2(v, u); % Angle in normalized coordinates
cos_theta = cos(theta);
sin_theta = sin(theta);

% Local radius of curvature approximation
% For ellipse, the local characteristic length varies with position
local_scale = sqrt((a * cos_theta).^2+(b * sin_theta).^2);

% More accurate distance using local scaling
% This provides better approximation near the ellipse boundary
d = (ellipse_val - 1) .* local_scale;

% Handle special case at origin to avoid division by zero
origin_idx = (abs(x) < eps) & (abs(y) < eps);
if any(origin_idx)
    d(origin_idx) = -min(a, b); % Distance from center to nearest boundary
end

end

