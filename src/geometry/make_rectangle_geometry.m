function G = make_rectangle_geometry(cfg)
%MAKE_RECTANGLE_GEOMETRY Generate mesh and boundaries for rectangle flow geometry
%
% This function encapsulates all rectangle-specific mesh generation and boundary
% extraction logic. It generates a rectangular domain with a rectangular
% obstacle using DistMesh, then extracts and classifies boundary nodes.
%
% INPUTS:
%   cfg - Configuration structure containing domain, mesh, and geometry parameters
%
% OUTPUTS:
%   G - Geometry structure containing:
%       .xy, .xy_s, .xt - Velocity grid nodes, pressure grid nodes, triangle connectivity
%       .boundary_in, .boundary_out, .boundary_y, .boundary_obs - Velocity boundary nodes
%       .boundary_in_s, .boundary_out_s, .boundary_y_s, .boundary_obs_s - Pressure boundary nodes
%       .rect_width, .rect_height - Rectangle dimensions
%       .rect_x_center, .rect_y_center - Rectangle center coordinates
%       .obs_normals_s - Unit normal vectors at obstacle boundary nodes

%% Extract domain and mesh parameters
x_min = cfg.domain.x_min;
x_max = cfg.domain.x_max;
y_min = cfg.domain.y_min;
y_max = cfg.domain.y_max;

dist = cfg.mesh.dist;
rect_width = cfg.geometry.rect_width;
rect_height = cfg.geometry.rect_height;
rect_x_center = cfg.geometry.rect_x_center;
rect_y_center = cfg.geometry.rect_y_center;
eps = cfg.mesh.boundary_eps;
a1 = cfg.mesh.refine_a1;
b1 = cfg.mesh.refine_b1;
a2 = cfg.mesh.refine_a2;
b2 = cfg.mesh.refine_b2;

%% Define rectangle-specific distance functions for DistMesh
% Rectangle bounds for the obstacle
rect_x1 = rect_x_center - rect_width/2;
rect_x2 = rect_x_center + rect_width/2;
rect_y1 = rect_y_center - rect_height/2;
rect_y2 = rect_y_center + rect_height/2;

% Use rounded rectangle with very small corner radius for near-rectangular shape
% Small radius prevents sharp corner convergence issues
corner_radius = min(rect_width, rect_height) * 0.03; % 3% rounding - minimal but smooth

fprintf('Using rounded rectangle (width=%.2f, height=%.2f, radius=%.3f)...\n', ...
        rect_width, rect_height, corner_radius);

% fd: Signed distance function for domain geometry (outer rectangle minus rounded rectangle)
fd = @(p) ddiff(drectangle(p, x_min, x_max, y_min, y_max), ...
               drounded_rectangle(p, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius));

% fd1: Edge length function for local mesh refinement
% Provides smaller elements near rectangle and in wake region behind it
fd1 = @(p) min(a1+b1*abs(drounded_rectangle(p, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius)), ...
               a2+b2*abs(dpoly(p, [rect_x2, rect_y_center; x_max, rect_y_center])));

% Define fixed corner points to maintain rectangular domain shape
fix = [x_min, y_min; x_min, y_max; x_max, y_max; x_max, y_min];

%% Generate triangular mesh using DistMesh algorithm
% xy_s: nodes for pressure (P-grid), xt: triangle connectivity
fprintf('Generating rectangle mesh with iteration limits...\n');
% Use robust DistMesh with iteration limits to prevent hanging
[xy_s, xt] = distmesh2d_robust(fd, fd1, dist, [x_min, y_min; x_max, y_max], fix);
fprintf('Mesh generation completed.\n');

%% Generate velocity grid nodes (V-grid) by adding edge midpoints
% This creates a staggered grid arrangement for better pressure-velocity coupling
xy = zeros(length(xt)*cfg.mesh.edge_multiplier, 2);
for j = 1:length(xt)
    % Add midpoint of each triangle edge to create velocity nodes
    xy((j - 1)*cfg.mesh.edge_multiplier+1, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 2), :)) / 2; % Edge 1-2
    xy((j - 1)*cfg.mesh.edge_multiplier+2, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 3), :)) / 2; % Edge 1-3
    xy((j - 1)*cfg.mesh.edge_multiplier+3, :) = (xy_s(xt(j, 2), :) + xy_s(xt(j, 3), :)) / 2; % Edge 2-3
end
xy = unique(xy, 'rows'); % Remove duplicate nodes

%% Extract and classify boundary nodes for pressure grid (P-grid)
% Remove corner nodes first to avoid duplicate boundary classification
idx_corners = (xy_s(:, 1) < x_min + eps & xy_s(:, 2) < y_min + eps) | ...
    (xy_s(:, 1) < x_min + eps & xy_s(:, 2) > y_max - eps) | ...
    (xy_s(:, 1) > x_max - eps & xy_s(:, 2) > y_max - eps) | ...
    (xy_s(:, 1) > x_max - eps & xy_s(:, 2) < y_min + eps);
xy_s(idx_corners, :) = [];

% Extract inlet boundary nodes (left wall, x = x_min)
idx_b_in = xy_s(:, 1) < x_min + eps;
boundary_in_s = xy_s(idx_b_in, :);
xy_s(idx_b_in, :) = [];

% Extract wall boundary nodes (top and bottom walls, y = y_min/y_max)
idx_b_y = xy_s(:, 2) > y_max - eps | xy_s(:, 2) < y_min + eps;
boundary_y_s = xy_s(idx_b_y, :);
xy_s(idx_b_y, :) = [];

% Extract outlet boundary nodes (right wall, x = x_max)
idx_b_out = xy_s(:, 1) > x_max - eps;
boundary_out_s = xy_s(idx_b_out, :);
xy_s(idx_b_out, :) = [];

% Extract rounded rectangle boundary nodes (obstacle boundary) - ROUNDED RECTANGLE
% Use rounded rectangle distance function for boundary detection
rect_dist_s = abs(drounded_rectangle(xy_s, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius));
idx_b_rect = rect_dist_s < eps;
boundary_obs_s = xy_s(idx_b_rect, :);
xy_s(idx_b_rect, :) = [];

%% Extract and classify boundary nodes for velocity grid (V-grid)
% Extract inlet boundary nodes
idx_b_in = xy(:, 1) < x_min + eps;
boundary_in = xy(idx_b_in, :);
xy(idx_b_in, :) = [];

% Extract wall boundary nodes
idx_b_y = xy(:, 2) > y_max - eps | xy(:, 2) < y_min + eps;
boundary_y = xy(idx_b_y, :);
xy(idx_b_y, :) = [];

% Extract outlet boundary nodes
idx_b_out = xy(:, 1) > x_max - eps;
boundary_out = xy(idx_b_out, :);
xy(idx_b_out, :) = [];

% Extract rounded rectangle boundary nodes (velocity grid) - ROUNDED RECTANGLE
% Use rounded rectangle distance function for boundary detection
rect_dist_v = abs(drounded_rectangle(xy, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius));
idx_b_rect = rect_dist_v < eps;
boundary_obs = xy(idx_b_rect, :);
xy(idx_b_rect, :) = [];

%% Compute unit normal vectors at obstacle boundary nodes (pressure grid)
% For rounded rectangle, compute gradient numerically for smooth normals
obs_normals_s = zeros(size(boundary_obs_s));

for i = 1:size(boundary_obs_s, 1)
    x_pt = boundary_obs_s(i, 1);
    y_pt = boundary_obs_s(i, 2);

    % Compute gradient numerically using finite differences
    delta = eps * 100;
    dx = (drounded_rectangle([x_pt + delta, y_pt], rect_x_center, rect_y_center, rect_width, rect_height, corner_radius) - ...
          drounded_rectangle([x_pt - delta, y_pt], rect_x_center, rect_y_center, rect_width, rect_height, corner_radius)) / (2 * delta);
    dy = (drounded_rectangle([x_pt, y_pt + delta], rect_x_center, rect_y_center, rect_width, rect_height, corner_radius) - ...
          drounded_rectangle([x_pt, y_pt - delta], rect_x_center, rect_y_center, rect_width, rect_height, corner_radius)) / (2 * delta);

    % Normalize to unit vector
    grad_norm = sqrt(dx^2 + dy^2);
    if grad_norm > eps
        obs_normals_s(i, :) = [dx, dy] / grad_norm;
    else
        obs_normals_s(i, :) = [1, 0]; % Default for degenerate case
    end
end

%% Precompute proximity indices for special RBF treatment
% Distance threshold for nodes requiring special treatment
r_dist = cfg.distances.x_min; % Reuse existing distance parameter

% Compute distance from rounded rectangle boundary for all nodes
% For rounded rectangle, use consistent distance function
rect_dist_V = abs(drounded_rectangle(xy, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius));
rect_dist_P = abs(drounded_rectangle(xy_s, rect_x_center, rect_y_center, rect_width, rect_height, corner_radius));

% Velocity grid: nodes near obstacle (within r_dist of surface)
idx_near_obs_V = rect_dist_V < r_dist;

% Pressure grid: nodes near obstacle (within r_dist of surface)
idx_near_obs_P = rect_dist_P < r_dist;

% Also precompute indices far from obstacle for high-order treatment
x_min_dist = cfg.distances.x_min;
x_max_dist = cfg.distances.x_max;
y_min_dist = cfg.distances.y_min;
y_max_dist = cfg.distances.y_max;

% Velocity grid: nodes far from all boundaries (high-order region)
idx_far_boundaries_V = rect_dist_V > r_dist & ...
    xy(:, 1) > x_min + x_min_dist & ...
    xy(:, 2) < y_max - y_max_dist & ...
    xy(:, 2) > y_min + y_min_dist & ...
    xy(:, 1) < x_max - x_max_dist;

% Pressure grid: nodes far from all boundaries (high-order region)
idx_far_boundaries_P = rect_dist_P > r_dist & ...
    xy_s(:, 1) > x_min + x_min_dist & ...
    xy_s(:, 2) < y_max - y_max_dist & ...
    xy_s(:, 2) > y_min + y_min_dist & ...
    xy_s(:, 1) < x_max - x_max_dist;

%% Assemble geometry structure
G.xy = xy; % Interior velocity nodes
G.xy_s = xy_s; % Interior pressure nodes
G.xt = xt; % Triangle connectivity

% Velocity grid boundaries
G.boundary_in = boundary_in; % Inlet boundary nodes (V-grid)
G.boundary_out = boundary_out; % Outlet boundary nodes (V-grid)
G.boundary_y = boundary_y; % Wall boundary nodes (V-grid)
G.boundary_obs = boundary_obs; % Obstacle boundary nodes (V-grid)

% Pressure grid boundaries
G.boundary_in_s = boundary_in_s; % Inlet boundary nodes (P-grid)
G.boundary_out_s = boundary_out_s; % Outlet boundary nodes (P-grid)
G.boundary_y_s = boundary_y_s; % Wall boundary nodes (P-grid)
G.boundary_obs_s = boundary_obs_s; % Obstacle boundary nodes (P-grid)

% Geometry-specific data (rounded rectangle)
G.rect_width = rect_width; % Rectangle width
G.rect_height = rect_height; % Rectangle height
G.rect_x_center = rect_x_center; % Rectangle center X
G.rect_y_center = rect_y_center; % Rectangle center Y
G.corner_radius = corner_radius; % Corner rounding radius used
G.obs_normals_s = obs_normals_s; % Unit normal vectors at obstacle boundary

% Precomputed proximity indices for special RBF treatment
G.idx_near_obs_V = idx_near_obs_V; % Velocity nodes near obstacle
G.idx_near_obs_P = idx_near_obs_P; % Pressure nodes near obstacle
G.idx_far_boundaries_V = idx_far_boundaries_V; % Velocity nodes far from boundaries
G.idx_far_boundaries_P = idx_far_boundaries_P; % Pressure nodes far from boundaries

end
