function G = make_airfoil_geometry(cfg)
    %MAKE_AIRFOIL_GEOMETRY Generate mesh and boundaries for NACA airfoil flow geometry
    %
    % This function encapsulates all airfoil-specific mesh generation and boundary
    % extraction logic. It generates a rectangular domain with a NACA 4-digit
    % airfoil obstacle using DistMesh, then extracts and classifies boundary nodes.
    %
    % INPUTS:
    %   cfg - Configuration structure containing domain, mesh, and geometry parameters
    %
    % OUTPUTS:
    %   G - Geometry structure containing:
    %       .xy, .xy_s, .xt - Velocity grid nodes, pressure grid nodes, triangle connectivity
    %       .boundary_in, .boundary_out, .boundary_y, .boundary_obs - Velocity boundary nodes
    %       .boundary_in_s, .boundary_out_s, .boundary_y_s, .boundary_obs_s - Pressure boundary nodes
    %       .naca_digits, .chord_length, .angle_of_attack - Airfoil parameters
    %       .airfoil_x_center, .airfoil_y_center - Airfoil center coordinates
    %       .obs_normals_s - Unit normal vectors at obstacle boundary nodes

    %% Extract domain and mesh parameters
    x_min = cfg.domain.x_min;
    x_max = cfg.domain.x_max;
    y_min = cfg.domain.y_min;
    y_max = cfg.domain.y_max;

    dist = cfg.mesh.dist;
    naca_digits = cfg.geometry.naca_digits;
    chord_length = cfg.geometry.chord_length;
    angle_of_attack = cfg.geometry.angle_of_attack;
    airfoil_x_center = cfg.geometry.airfoil_x_center;
    airfoil_y_center = cfg.geometry.airfoil_y_center;
    eps = cfg.mesh.boundary_eps;
    a1 = cfg.mesh.refine_a1;
    b1 = cfg.mesh.refine_b1;
    a2 = cfg.mesh.refine_a2;
    b2 = cfg.mesh.refine_b2;

    %% Define airfoil-specific distance functions for DistMesh
    % fd: Signed distance function for domain geometry (rectangle minus airfoil)
    fd = @(p) ddiff(drectangle(p, x_min, x_max, y_min, y_max), ...
                    dairfoil(p, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack));

    % fd1: Edge length function for local mesh refinement
    % Provides smaller elements near airfoil and in wake region behind it
    % For wake refinement, use a line from trailing edge to domain outlet
    trailing_edge_x = airfoil_x_center + chord_length * cos(angle_of_attack * pi / 180);
    trailing_edge_y = airfoil_y_center + chord_length * sin(angle_of_attack * pi / 180);

    fd1 = @(p) min(a1 + b1 * abs(dairfoil(p, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack)), ...
                   a2 + b2 * abs(dpoly(p, [trailing_edge_x, trailing_edge_y; x_max, trailing_edge_y])));

    % Define fixed corner points to maintain rectangular domain shape
    fix = [x_min, y_min; x_min, y_max; x_max, y_max; x_max, y_min];

    %% Generate triangular mesh using DistMesh algorithm
    % xy_s: nodes for pressure (P-grid), xt: triangle connectivity
    fprintf('Generating NACA %d%d%d%d airfoil mesh (chord=%.2f, AoA=%.1f deg)...\n', ...
            naca_digits(1), naca_digits(2), naca_digits(3), naca_digits(4), ...
            chord_length, angle_of_attack);

    % Use standard DistMesh for airfoil geometry
    [xy_s, xt] = distmesh2d(fd, fd1, dist, [x_min, y_min; x_max, y_max], fix);
    fprintf('Airfoil mesh generation completed.\n');

    %% Generate velocity grid nodes (V-grid) by adding edge midpoints
    % This creates a staggered grid arrangement for better pressure-velocity coupling
    % OPTIMIZED: Vectorized computation instead of loop
    
    n_triangles = size(xt, 1);
    edge_mult = cfg.mesh.edge_multiplier;
    
    % Pre-allocate for all edge midpoints
    xy = zeros(n_triangles * edge_mult, 2);
    
    % Vectorized computation of all edge midpoints
    % Edge 1-2 midpoints
    idx1 = 1:edge_mult:n_triangles * edge_mult;
    xy(idx1, :) = (xy_s(xt(:, 1), :) + xy_s(xt(:, 2), :)) / 2;
    
    % Edge 1-3 midpoints  
    idx2 = 2:edge_mult:n_triangles * edge_mult;
    xy(idx2, :) = (xy_s(xt(:, 1), :) + xy_s(xt(:, 3), :)) / 2;
    
    % Edge 2-3 midpoints
    idx3 = 3:edge_mult:n_triangles * edge_mult;
    xy(idx3, :) = (xy_s(xt(:, 2), :) + xy_s(xt(:, 3), :)) / 2;
    
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

    % Extract airfoil boundary nodes (airfoil boundary) - AIRFOIL-SPECIFIC
    % Use airfoil distance function for boundary detection
    airfoil_dist_s = abs(dairfoil(xy_s, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack));
    idx_b_airfoil = airfoil_dist_s < eps;
    boundary_obs_s = xy_s(idx_b_airfoil, :);
    xy_s(idx_b_airfoil, :) = [];

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

    % Extract airfoil boundary nodes - AIRFOIL-SPECIFIC
    % Use same airfoil distance function for velocity grid
    airfoil_dist_v_full = abs(dairfoil(xy, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack));
    idx_b_airfoil = airfoil_dist_v_full < eps;
    boundary_obs = xy(idx_b_airfoil, :);
    xy(idx_b_airfoil, :) = [];

    %% Compute obstacle normal vectors for pressure boundary conditions
    % For airfoil: compute normal vectors numerically using gradient of distance function
    % OPTIMIZED: Vectorized computation to avoid expensive loop
    
    n_boundary = size(boundary_obs_s, 1);
    obs_normals_s = zeros(n_boundary, 2);
    
    if n_boundary > 0
        % Compute gradient numerically using finite differences - VECTORIZED
        delta = eps * 100;
        
        % Create perturbed point arrays for vectorized dairfoil calls
        x_pts = boundary_obs_s(:, 1);
        y_pts = boundary_obs_s(:, 2);
        
        % Compute x-derivatives (4 calls instead of 4*N calls)
        pts_x_plus = [x_pts + delta, y_pts];
        pts_x_minus = [x_pts - delta, y_pts];
        pts_y_plus = [x_pts, y_pts + delta];
        pts_y_minus = [x_pts, y_pts - delta];
        
        % Vectorized dairfoil calls
        d_x_plus = dairfoil(pts_x_plus, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack);
        d_x_minus = dairfoil(pts_x_minus, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack);
        d_y_plus = dairfoil(pts_y_plus, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack);
        d_y_minus = dairfoil(pts_y_minus, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack);
        
        % Compute gradients
        dx = (d_x_plus - d_x_minus) / (2 * delta);
        dy = (d_y_plus - d_y_minus) / (2 * delta);
        
        % Normalize to unit vectors - VECTORIZED
        grad_norm = sqrt(dx.^2 + dy.^2);
        valid_idx = grad_norm > eps;
        
        % Set valid normals
        obs_normals_s(valid_idx, 1) = dx(valid_idx) ./ grad_norm(valid_idx);
        obs_normals_s(valid_idx, 2) = dy(valid_idx) ./ grad_norm(valid_idx);
        
        % Set default for degenerate cases
        obs_normals_s(~valid_idx, :) = repmat([1, 0], sum(~valid_idx), 1);
    end

    % Verify normal vectors are reasonable
    if any(isnan(obs_normals_s(:))) || any(isinf(obs_normals_s(:)))
        error('Invalid normal vectors computed for airfoil boundary');
    end

    %% Precompute near-obstacle node indices for special RBF treatment
    % Distance threshold for near-obstacle treatment
    r_dist = cfg.mesh.refine_a2 * cfg.mesh.edge_multiplier;

    % For airfoil, use distance to airfoil boundary for proximity calculations
    % OPTIMIZED: Compute distance for final xy array (after boundary extraction)
    airfoil_dist_V = abs(dairfoil(xy, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack));
    airfoil_dist_P = abs(dairfoil(xy_s, airfoil_x_center, airfoil_y_center, naca_digits, chord_length, angle_of_attack));

    % Velocity grid: nodes near obstacle (within r_dist of surface)
    idx_near_obs_V = airfoil_dist_V < r_dist;

    % Pressure grid: nodes near obstacle (within r_dist of surface)
    idx_near_obs_P = airfoil_dist_P < r_dist;

    % Also precompute indices far from obstacle for high-order treatment
    x_min_dist = cfg.distances.x_min;
    x_max_dist = cfg.distances.x_max;
    y_min_dist = cfg.distances.y_min;
    y_max_dist = cfg.distances.y_max;

    % Velocity grid: nodes far from all boundaries (high-order region)
    idx_far_boundaries_V = airfoil_dist_V > r_dist & ...
      xy(:, 1) > x_min + x_min_dist & ...
      xy(:, 2) < y_max - y_max_dist & ...
      xy(:, 2) > y_min + y_min_dist & ...
      xy(:, 1) < x_max - x_max_dist;

    % Pressure grid: nodes far from all boundaries (high-order region)
    idx_far_boundaries_P = airfoil_dist_P > r_dist & ...
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
    G.naca_digits = naca_digits; % NACA 4-digit parameters
    G.chord_length = chord_length; % Airfoil chord length
    G.angle_of_attack = angle_of_attack; % Angle of attack in degrees
    G.airfoil_x_center = airfoil_x_center; % Airfoil center X
    G.airfoil_y_center = airfoil_y_center; % Airfoil center Y
    G.obs_normals_s = obs_normals_s; % Unit normal vectors at obstacle boundary

    % Precomputed proximity indices for special RBF treatment
    G.idx_near_obs_V = idx_near_obs_V; % Velocity nodes near obstacle
    G.idx_near_obs_P = idx_near_obs_P; % Pressure nodes near obstacle
    G.idx_far_boundaries_V = idx_far_boundaries_V; % Velocity nodes far from boundaries
    G.idx_far_boundaries_P = idx_far_boundaries_P; % Pressure nodes far from boundaries

end
