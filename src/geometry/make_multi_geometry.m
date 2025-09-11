function G = make_multi_geometry(cfg)
  %MAKE_MULTI_GEOMETRY Generate mesh and boundaries for multiple obstacle flow geometry
  %
  % This function encapsulates all multi-obstacle mesh generation and boundary
  % extraction logic. It generates a rectangular domain with multiple obstacles
  % using DistMesh, then extracts and classifies boundary nodes.
  %
  % INPUTS:
  %   cfg - Configuration structure containing domain, mesh, and geometry parameters
  %         cfg.geometry.obstacles - Array of obstacle structures, each containing:
  %           .type   - Obstacle type ('cylinder', 'ellipse', 'rectangle', 'airfoil')
  %           .center - [x, y] center coordinates
  %           .params - Structure with type-specific parameters
  %
  % OUTPUTS:
  %   G - Geometry structure containing:
  %       .xy, .xy_s, .xt - Velocity grid nodes, pressure grid nodes, triangle connectivity
  %       .boundary_in, .boundary_out, .boundary_y, .boundary_obs - Velocity boundary nodes
  %       .boundary_in_s, .boundary_out_s, .boundary_y_s, .boundary_obs_s - Pressure boundary nodes
  %       .obs_normals_s - Unit normal vectors at obstacle boundary nodes
  %       .boundary_obs_ids - Obstacle ID for each boundary node (for debugging)

  %% Extract domain and mesh parameters
  x_min = cfg.domain.x_min;
  x_max = cfg.domain.x_max;
  y_min = cfg.domain.y_min;
  y_max = cfg.domain.y_max;

  dist = cfg.mesh.dist;
  eps = cfg.mesh.boundary_eps;
  a1 = cfg.mesh.refine_a1;
  b1 = cfg.mesh.refine_b1;
  a2 = cfg.mesh.refine_a2;
  b2 = cfg.mesh.refine_b2;

  obstacles = cfg.geometry.obstacles;
  num_obstacles = length(obstacles);

  fprintf('Building multi-obstacle geometry with %d obstacles:\n', num_obstacles);
  for i = 1:num_obstacles
    obs = obstacles(i);
    fprintf('  %d: %s at [%.2f, %.2f]\n', i, obs.type, obs.center(1), obs.center(2));
  end

  %% Build union signed distance function for all obstacles
  fd_obs_union = build_obstacles_union_sdf(obstacles);

  %% Define domain and mesh refinement functions
  % Domain: rectangular domain minus union of obstacles
  fd = @(p) ddiff(drectangle(p, x_min, x_max, y_min, y_max), fd_obs_union(p));

  % Mesh refinement: smaller elements near obstacles and individual wake regions
  % Create wake refinement for each obstacle individually
  wake_refinements = cell(length(obstacles), 1);
  for i = 1:length(obstacles)
    obs = obstacles(i);
    center = obs.center;
    % Find the downstream edge of each obstacle
    switch lower(obs.type)
      case 'cylinder'
        wake_start_x = center(1) + obs.params.radius;
        wake_start_y = center(2);
      case 'ellipse'
        wake_start_x = center(1) + obs.params.a;
        wake_start_y = center(2);
      case 'rectangle'
        wake_start_x = center(1) + obs.params.width / 2;
        wake_start_y = center(2);
      otherwise
        wake_start_x = center(1) + 0.5;
        wake_start_y = center(2);
    end
    % Create wake line for this obstacle
    wake_refinements{i} = @(p) a2 + b2 * abs(dpoly(p, [wake_start_x, wake_start_y; x_max, wake_start_y]));
  end

  % Combine all refinements: obstacle + all individual wakes
  fd1 = @(p) compute_combined_refinement(p, fd_obs_union, wake_refinements, a1, b1);

  % Define fixed corner points to maintain rectangular domain shape
  fix = [x_min, y_min; x_min, y_max; x_max, y_max; x_max, y_min];

  %% Generate triangular mesh using robust DistMesh algorithm
  fprintf('Generating multi-obstacle mesh with iteration limits...\n');
  [xy_s, xt] = distmesh2d_robust(fd, fd1, dist, [x_min, y_min; x_max, y_max], fix);
  fprintf('Multi-obstacle mesh generation completed.\n');

  %% Generate velocity grid nodes (V-grid) by adding edge midpoints
  xy = compute_velocity_grid_from_triangles(xy_s, xt, cfg.mesh.edge_multiplier);

  %% Extract and classify boundary nodes for pressure grid (P-grid)
  [boundary_in_s, boundary_y_s, boundary_out_s, boundary_obs_s, xy_s, obs_ids_s] = ...
      classify_pressure_boundaries(xy_s, x_min, x_max, y_min, y_max, fd_obs_union, obstacles, eps);

  %% Extract and classify boundary nodes for velocity grid (V-grid)
  % Use larger tolerance for velocity grid obstacle boundaries since edge midpoints
  % are typically further from obstacle surfaces than pressure grid nodes
  eps_v = eps * 2.5; % Increased tolerance for velocity grid obstacle boundaries
  [boundary_in, boundary_y, boundary_out, boundary_obs, xy, ~] = ...
      classify_velocity_boundaries(xy, x_min, x_max, y_min, y_max, fd_obs_union, obstacles, eps_v);

  %% Compute unit normal vectors at obstacle boundary nodes (pressure grid)
  obs_normals_s = compute_obstacle_normals(fd_obs_union, boundary_obs_s, eps);

  %% Precompute proximity indices for special RBF treatment
  r_dist = cfg.distances.x_min; % Distance threshold for near-obstacle treatment

  % Compute distance from obstacle union for all nodes
  obs_dist_V = abs(fd_obs_union(xy));
  obs_dist_P = abs(fd_obs_union(xy_s));

  % Velocity grid: nodes near obstacles (within r_dist of surface)
  idx_near_obs_V = obs_dist_V < r_dist;

  % Pressure grid: nodes near obstacles (within r_dist of surface)
  idx_near_obs_P = obs_dist_P < r_dist;

  % Distance thresholds for far-boundary treatment
  x_min_dist = cfg.distances.x_min;
  x_max_dist = cfg.distances.x_max;
  y_min_dist = cfg.distances.y_min;
  y_max_dist = cfg.distances.y_max;

  % Velocity grid: nodes far from all boundaries (high-order region)
  idx_far_boundaries_V = obs_dist_V > r_dist & ...
      xy(:, 1) > x_min + x_min_dist & ...
      xy(:, 2) < y_max - y_max_dist & ...
      xy(:, 2) > y_min + y_min_dist & ...
      xy(:, 1) < x_max - x_max_dist;

  % Pressure grid: nodes far from all boundaries (high-order region)
  idx_far_boundaries_P = obs_dist_P > r_dist & ...
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

  % Obstacle-specific data
  G.obs_normals_s = obs_normals_s; % Unit normal vectors at obstacle boundary
  G.boundary_obs_ids = obs_ids_s; % Obstacle ID for each boundary node (debugging)

  % Precomputed proximity indices for special RBF treatment
  G.idx_near_obs_V = idx_near_obs_V; % Velocity nodes near obstacles
  G.idx_near_obs_P = idx_near_obs_P; % Pressure nodes near obstacles
  G.idx_far_boundaries_V = idx_far_boundaries_V; % Velocity nodes far from boundaries
  G.idx_far_boundaries_P = idx_far_boundaries_P; % Pressure nodes far from boundaries

  % Store obstacle configuration for reference
  G.obstacles = obstacles;
  G.num_obstacles = num_obstacles;

  fprintf('Multi-obstacle geometry completed: %d obstacles, %d P-nodes, %d V-nodes\n', ...
          num_obstacles, length(xy_s), length(xy));
end

function fd_obs_union = build_obstacles_union_sdf(obstacles)
  %BUILD_OBSTACLES_UNION_SDF Build union signed distance function for all obstacles
  %
  % INPUTS:
  %   obstacles - Array of obstacle structures
  %
  % OUTPUTS:
  %   fd_obs_union - Function handle for union SDF

  num_obstacles = length(obstacles);
  if num_obstacles == 0
    error('At least one obstacle must be specified for multi geometry');
  end

  % Build individual SDFs
  obstacle_sdfs = cell(num_obstacles, 1);
  for i = 1:num_obstacles
    obs = obstacles(i);
    type = lower(obs.type);
    center = obs.center;
    params = obs.params;

    switch type
      case 'cylinder'
        obstacle_sdfs{i} = @(p) dcircle(pshift(p, -center(1), -center(2)), 0, 0, params.radius);

      case 'ellipse'
        if isfield(params, 'theta_deg') && params.theta_deg ~= 0
          % Rotated ellipse
          theta_rad = params.theta_deg * pi / 180;
          obstacle_sdfs{i} = @(p) dellipse(protate(pshift(p, -center(1), -center(2)), -theta_rad), ...
                                           0, 0, params.a, params.b);
        else
          % Axis-aligned ellipse
          obstacle_sdfs{i} = @(p) dellipse(pshift(p, -center(1), -center(2)), 0, 0, params.a, params.b);
        end

      case 'rectangle'
        % Use rounded rectangle with small corner radius
        corner_radius_frac = 0.03; % Default
        if isfield(params, 'corner_radius_frac')
          corner_radius_frac = params.corner_radius_frac;
        end
        corner_radius = min(params.width, params.height) * corner_radius_frac;
        obstacle_sdfs{i} = @(p) drounded_rectangle(p, center(1), center(2), ...
                                                   params.width, params.height, corner_radius);

      case 'airfoil'
        obstacle_sdfs{i} = @(p) dairfoil(p, center(1), center(2), params.naca_digits, ...
                                         params.chord_length, params.angle_of_attack);

      otherwise
        error('Unsupported obstacle type: %s', type);
    end
  end

  % Build union using recursive dunion
  fd_obs_union = obstacle_sdfs{1};
  for i = 2:num_obstacles
    current_sdf = obstacle_sdfs{i};
    prev_sdf = fd_obs_union; % Capture current union for closure
    fd_obs_union = @(p) min(prev_sdf(p), current_sdf(p));
  end
end

function xy = compute_velocity_grid_from_triangles(xy_s, xt, edge_multiplier)
  %COMPUTE_VELOCITY_GRID_FROM_TRIANGLES Generate velocity grid from triangle edges
  %
  % INPUTS:
  %   xy_s - Pressure grid nodes
  %   xt   - Triangle connectivity
  %   edge_multiplier - Number of edges per triangle (should be 3)
  %
  % OUTPUTS:
  %   xy - Velocity grid nodes (edge midpoints)

  % Generate velocity grid nodes (V-grid) by adding edge midpoints
  % This creates a staggered grid arrangement for better pressure-velocity coupling
  xy = zeros(length(xt) * edge_multiplier, 2);
  for j = 1:length(xt)
    % Add midpoint of each triangle edge to create velocity nodes
    xy((j - 1) * edge_multiplier + 1, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 2), :)) / 2; % Edge 1-2
    xy((j - 1) * edge_multiplier + 2, :) = (xy_s(xt(j, 1), :) + xy_s(xt(j, 3), :)) / 2; % Edge 1-3
    xy((j - 1) * edge_multiplier + 3, :) = (xy_s(xt(j, 2), :) + xy_s(xt(j, 3), :)) / 2; % Edge 2-3
  end
  xy = unique(xy, 'rows'); % Remove duplicate nodes
end

function [boundary_in_s, boundary_y_s, boundary_out_s, boundary_obs_s, xy_s, obs_ids_s] = ...
    classify_pressure_boundaries(xy_s, x_min, x_max, y_min, y_max, fd_obs_union, obstacles, eps)
  %CLASSIFY_PRESSURE_BOUNDARIES Extract and classify boundary nodes for pressure grid
  %
  % INPUTS:
  %   xy_s - Pressure grid nodes
  %   x_min, x_max, y_min, y_max - Domain boundaries
  %   fd_obs_union - Union SDF for obstacles
  %   obstacles - Obstacle array (for ID assignment)
  %   eps - Boundary tolerance
  %
  % OUTPUTS:
  %   boundary_*_s - Boundary node arrays
  %   xy_s - Interior nodes (boundaries removed)
  %   obs_ids_s - Obstacle IDs for obstacle boundary nodes

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

  % Extract obstacle boundary nodes using union SDF
  obs_dist_s = abs(fd_obs_union(xy_s));
  idx_b_obs = obs_dist_s < eps;
  boundary_obs_s = xy_s(idx_b_obs, :);
  xy_s(idx_b_obs, :) = [];

  % Assign obstacle IDs to boundary nodes (for debugging/visualization)
  obs_ids_s = assign_obstacle_ids(boundary_obs_s, obstacles);
end

function [boundary_in, boundary_y, boundary_out, boundary_obs, xy, obs_ids] = ...
    classify_velocity_boundaries(xy, x_min, x_max, y_min, y_max, fd_obs_union, obstacles, eps)
  %CLASSIFY_VELOCITY_BOUNDARIES Extract and classify boundary nodes for velocity grid
  %
  % INPUTS:
  %   xy - Velocity grid nodes
  %   x_min, x_max, y_min, y_max - Domain boundaries
  %   fd_obs_union - Union SDF for obstacles
  %   obstacles - Obstacle array (for ID assignment)
  %   eps - Boundary tolerance
  %
  % OUTPUTS:
  %   boundary_* - Boundary node arrays
  %   xy - Interior nodes (boundaries removed)
  %   obs_ids - Obstacle IDs for obstacle boundary nodes

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

  % Extract obstacle boundary nodes using union SDF
  obs_dist = abs(fd_obs_union(xy));
  idx_b_obs = obs_dist < eps;
  boundary_obs = xy(idx_b_obs, :);
  xy(idx_b_obs, :) = [];

  % Assign obstacle IDs to boundary nodes (for debugging/visualization)
  obs_ids = assign_obstacle_ids(boundary_obs, obstacles);
end

function obs_ids = assign_obstacle_ids(boundary_pts, obstacles)
  %ASSIGN_OBSTACLE_IDS Assign obstacle ID to each boundary point
  %
  % INPUTS:
  %   boundary_pts - Boundary node coordinates [N x 2]
  %   obstacles - Obstacle array
  %
  % OUTPUTS:
  %   obs_ids - Obstacle ID for each boundary point [N x 1]

  num_pts = size(boundary_pts, 1);
  num_obstacles = length(obstacles);
  obs_ids = zeros(num_pts, 1);

  if num_pts == 0
    return
  end

  % For each boundary point, find the closest obstacle
  for i = 1:num_pts
    pt = boundary_pts(i, :);
    min_dist = inf;
    closest_obs = 1;

    for j = 1:num_obstacles
      obs = obstacles(j);
      % Compute distance to this obstacle
      switch lower(obs.type)
        case 'cylinder'
          dist = abs(dcircle(pshift(pt, -obs.center(1), -obs.center(2)), 0, 0, obs.params.radius));
        case 'ellipse'
          if isfield(obs.params, 'theta_deg') && obs.params.theta_deg ~= 0
            theta_rad = obs.params.theta_deg * pi / 180;
            dist = abs(dellipse(protate(pshift(pt, -obs.center(1), -obs.center(2)), -theta_rad), ...
                                0, 0, obs.params.a, obs.params.b));
          else
            dist = abs(dellipse(pshift(pt, -obs.center(1), -obs.center(2)), 0, 0, obs.params.a, obs.params.b));
          end
        case 'rectangle'
          corner_radius_frac = 0.03;
          if isfield(obs.params, 'corner_radius_frac')
            corner_radius_frac = obs.params.corner_radius_frac;
          end
          corner_radius = min(obs.params.width, obs.params.height) * corner_radius_frac;
          dist = abs(drounded_rectangle(pt, obs.center(1), obs.center(2), ...
                                        obs.params.width, obs.params.height, corner_radius));
        case 'airfoil'
          dist = abs(dairfoil(pt, obs.center(1), obs.center(2), obs.params.naca_digits, ...
                              obs.params.chord_length, obs.params.angle_of_attack));
        otherwise
          dist = inf;
      end

      if dist < min_dist
        min_dist = dist;
        closest_obs = j;
      end
    end

    obs_ids(i) = closest_obs;
  end
end

function normals = compute_obstacle_normals(fd_obs_union, boundary_pts, eps)
  %COMPUTE_OBSTACLE_NORMALS Compute unit normal vectors at obstacle boundary points
  %
  % INPUTS:
  %   fd_obs_union - Union SDF for obstacles
  %   boundary_pts - Boundary point coordinates [N x 2]
  %   eps - Numerical differentiation step size
  %
  % OUTPUTS:
  %   normals - Unit normal vectors [N x 2]

  num_pts = size(boundary_pts, 1);
  normals = zeros(num_pts, 2);

  if num_pts == 0
    return
  end

  % Compute numerical gradient using finite differences (vectorized)
  delta = eps * 100;

  % Create perturbed point arrays for vectorized SDF calls
  x_pts = boundary_pts(:, 1);
  y_pts = boundary_pts(:, 2);

  % Compute finite difference gradients
  pts_x_plus = [x_pts + delta, y_pts];
  pts_x_minus = [x_pts - delta, y_pts];
  pts_y_plus = [x_pts, y_pts + delta];
  pts_y_minus = [x_pts, y_pts - delta];

  % Vectorized SDF calls
  d_x_plus = fd_obs_union(pts_x_plus);
  d_x_minus = fd_obs_union(pts_x_minus);
  d_y_plus = fd_obs_union(pts_y_plus);
  d_y_minus = fd_obs_union(pts_y_minus);

  % Compute gradients
  dx = (d_x_plus - d_x_minus) / (2 * delta);
  dy = (d_y_plus - d_y_minus) / (2 * delta);

  % Normalize to unit vectors
  grad_norm = sqrt(dx.^2 + dy.^2);
  valid_idx = grad_norm > eps;

  % Set valid normals
  normals(valid_idx, 1) = dx(valid_idx) ./ grad_norm(valid_idx);
  normals(valid_idx, 2) = dy(valid_idx) ./ grad_norm(valid_idx);

  % Set default for degenerate cases
  normals(~valid_idx, :) = repmat([1, 0], sum(~valid_idx), 1);

  % Verify normal vectors are reasonable
  if any(isnan(normals(:))) || any(isinf(normals(:)))
    error('Invalid normal vectors computed for multi-obstacle boundary');
  end
end

function refinement = compute_combined_refinement(p, fd_obs_union, wake_refinements, a1, b1)
  %COMPUTE_COMBINED_REFINEMENT Combine obstacle and individual wake refinements
  %
  % INPUTS:
  %   p                - Points to evaluate [N x 2]
  %   fd_obs_union     - Union of all obstacle distance functions
  %   wake_refinements - Cell array of wake refinement functions
  %   a1, b1           - Obstacle refinement parameters
  %
  % OUTPUTS:
  %   refinement       - Combined refinement values [N x 1]

  % Start with obstacle refinement
  obstacle_ref = a1 + b1 * abs(fd_obs_union(p));

  % Evaluate all wake refinements
  num_wakes = length(wake_refinements);
  if num_wakes == 0
    refinement = obstacle_ref;
    return
  end

  % Compute all wake refinements and take minimum with obstacle refinement
  all_refinements = zeros(size(p, 1), num_wakes + 1);
  all_refinements(:, 1) = obstacle_ref;

  for i = 1:num_wakes
    all_refinements(:, i + 1) = wake_refinements{i}(p);
  end

  % Take minimum across all refinement types (obstacle + all wakes)
  refinement = min(all_refinements, [], 2);
end
