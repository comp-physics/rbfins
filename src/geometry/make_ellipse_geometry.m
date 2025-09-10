function G = make_ellipse_geometry(cfg)
  %MAKE_ELLIPSE_GEOMETRY Generate mesh and boundaries for ellipse flow geometry
  %
  % This function encapsulates all ellipse-specific mesh generation and boundary
  % extraction logic. It generates a rectangular domain with an elliptical
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
  %       .a, .b - Ellipse semi-major and semi-minor axes
  %       .obs_normals_s - Unit normal vectors at obstacle boundary nodes
  %% Extract domain and mesh parameters
  x_min = cfg.domain.x_min;
  x_max = cfg.domain.x_max;
  y_min = cfg.domain.y_min;
  y_max = cfg.domain.y_max;

  dist = cfg.mesh.dist;
  a = cfg.geometry.ellipse_a; % Semi-major axis (x-direction)
  b = cfg.geometry.ellipse_b; % Semi-minor axis (y-direction)
  eps = cfg.mesh.boundary_eps;
  a1 = cfg.mesh.refine_a1;
  b1 = cfg.mesh.refine_b1;
  a2 = cfg.mesh.refine_a2;
  b2 = cfg.mesh.refine_b2;
  %% Define ellipse-specific distance functions for DistMesh
  % fd: Signed distance function for domain geometry (rectangle minus ellipse)
  fd = @(p) ddiff(drectangle(p, x_min, x_max, y_min, y_max), dellipse(p, 0, 0, a, b));

  % fd1: Edge length function for local mesh refinement
  % Provides smaller elements near ellipse and in wake region behind it
  fd1 = @(p) min(a1 + b1 * abs(dellipse(p, 0, 0, a, b)), a2 + b2 * abs(dpoly(p, [a, 0; x_max, 0])));

  % Define fixed corner points to maintain rectangular domain shape
  fix = [x_min, y_min; x_min, y_max; x_max, y_max; x_max, y_min];
  %% Generate triangular mesh using DistMesh algorithm
  % xy_s: nodes for pressure (P-grid), xt: triangle connectivity
  [xy_s, xt] = distmesh2d(fd, fd1, dist, [x_min, y_min; x_max, y_max], fix);
  %% Generate velocity grid nodes (V-grid) by adding edge midpoints
  % This creates a staggered grid arrangement for better pressure-velocity coupling
  xy = zeros(length(xt) * cfg.mesh.edge_multiplier, 2);
  for j = 1:length(xt)
    % Add midpoint of each triangle edge to create velocity nodes
    xy((j - 1) * cfg.mesh.edge_multiplier + 1, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 2), :)) / 2; % Edge 1-2
    xy((j - 1) * cfg.mesh.edge_multiplier + 2, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 3), :)) / 2; % Edge 1-3
    xy((j - 1) * cfg.mesh.edge_multiplier + 3, :) = (xy_s(xt(j, 2), :) + xy_s(xt(j, 3), :)) / 2; % Edge 2-3
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

  % Extract ellipse boundary nodes (elliptical boundary) - ELLIPSE-SPECIFIC
  % Use more robust boundary detection with proper scaling
  ellipse_tol = eps * 5; % More conservative tolerance
  ellipse_dist_s = abs(sqrt((xy_s(:, 1) / a).^2 + (xy_s(:, 2) / b).^2) - 1);
  idx_b_e = ellipse_dist_s < ellipse_tol;
  boundary_obs_s = xy_s(idx_b_e, :);
  xy_s(idx_b_e, :) = [];
  %% Extract and classify boundary nodes for velocity grid (V-grid)
  % Extract inlet boundary nodes
  idx_b_in = xy(:, 1) < x_min + eps;
  boundary_in = xy(idx_b_in, :);
  xy(idx_b_in, :) = [];

  % Extract wall boundary nodes (top and bottom)
  idx_b_y = xy(:, 2) > y_max - eps | xy(:, 2) < y_min + eps;
  boundary_y = xy(idx_b_y, :);
  xy(idx_b_y, :) = [];

  % Extract outlet boundary nodes
  idx_b_out = xy(:, 1) > x_max - eps;
  boundary_out = xy(idx_b_out, :);
  xy(idx_b_out, :) = [];

  % Extract ellipse boundary nodes - ELLIPSE-SPECIFIC
  % Use same robust boundary detection for velocity grid
  ellipse_dist_v = abs(sqrt((xy(:, 1) / a).^2 + (xy(:, 2) / b).^2) - 1);
  idx_b_e = ellipse_dist_v < ellipse_tol;
  boundary_obs = xy(idx_b_e, :);
  xy(idx_b_e, :) = [];
  %% Compute obstacle normal vectors for pressure boundary conditions
  % For ellipse: normal vector at point (x,y) is proportional to [x/a^2, y/b^2]
  % First compute the unnormalized normals
  normal_x = boundary_obs_s(:, 1) / (a^2);
  normal_y = boundary_obs_s(:, 2) / (b^2);
  normal_mag = sqrt(normal_x.^2 + normal_y.^2);

  % Add small regularization to avoid division by zero
  normal_mag = max(normal_mag, eps);

  % Unit normals (outward pointing)
  obs_normals_s = [normal_x ./ normal_mag, normal_y ./ normal_mag];

  % Verify normal vectors are reasonable
  if any(isnan(obs_normals_s(:))) || any(isinf(obs_normals_s(:)))
    error('Invalid normal vectors computed for ellipse boundary');
  end
  %% Precompute near-obstacle node indices for special RBF treatment
  % Distance threshold for near-obstacle treatment
  r_dist = cfg.mesh.refine_a2 * cfg.mesh.edge_multiplier;

  % For ellipse, use distance to ellipse boundary for proximity calculations
  % Use improved distance calculation consistent with dellipse function
  ellipse_dist_V = abs(dellipse(xy, 0, 0, a, b));
  ellipse_dist_P = abs(dellipse(xy_s, 0, 0, a, b));

  % Velocity grid: nodes near obstacle (within r_dist of surface)
  idx_near_obs_V = ellipse_dist_V < r_dist;

  % Pressure grid: nodes near obstacle (within r_dist of surface)
  idx_near_obs_P = ellipse_dist_P < r_dist;

  % Also precompute indices far from obstacle for high-order treatment
  x_min_dist = cfg.distances.x_min;
  x_max_dist = cfg.distances.x_max;
  y_min_dist = cfg.distances.y_min;
  y_max_dist = cfg.distances.y_max;

  % Velocity grid: nodes far from all boundaries (high-order region)
  idx_far_boundaries_V = ellipse_dist_V > r_dist & ...
    xy(:, 1) > x_min + x_min_dist & ...
    xy(:, 2) < y_max - y_max_dist & ...
    xy(:, 2) > y_min + y_min_dist & ...
    xy(:, 1) < x_max - x_max_dist;

  % Pressure grid: nodes far from all boundaries (high-order region)
  idx_far_boundaries_P = ellipse_dist_P > r_dist & ...
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

  % Geometry-specific data
  G.a = a; % Ellipse semi-major axis
  G.b = b; % Ellipse semi-minor axis
  G.obs_normals_s = obs_normals_s; % Unit normal vectors at obstacle boundary

  % Precomputed proximity indices for special RBF treatment
  G.idx_near_obs_V = idx_near_obs_V; % Velocity nodes near obstacle
  G.idx_near_obs_P = idx_near_obs_P; % Pressure nodes near obstacle
  G.idx_far_boundaries_V = idx_far_boundaries_V; % Velocity nodes far from boundaries
  G.idx_far_boundaries_P = idx_far_boundaries_P; % Pressure nodes far from boundaries

end
