clc;
clear;

% Add required paths for src functions and lib dependencies
scriptDir = fileparts(mfilename('fullpath'));

% Add src directory for supporting functions (config.m is now in main directory)
srcPath = fullfile(scriptDir, 'src');
if exist(srcPath, 'dir')
  addpath(srcPath);
  % Also add geometry subdirectory
  geometryPath = fullfile(srcPath, 'geometry');
  if exist(geometryPath, 'dir')
    addpath(geometryPath);
  end
end

% Load configuration parameters
% Can be overridden by setting CONFIG_FUNC and GEOMETRY_TYPE environment variables
config_func_name = getenv('CONFIG_FUNC');
geometry_type = getenv('GEOMETRY_TYPE');

if isempty(config_func_name)
  config_func_name = 'config';
end

if isempty(geometry_type)
  % Default behavior - call config function without arguments
  cfg = config();
else
  % Pass geometry type to config function
  cfg = config(geometry_type);
end

% Make config globally accessible for diagnostic functions
global cfg

%% 1) Setup environment and dependencies
[doPlot, isCI, isTest, Nt] = setup_environment(cfg, scriptDir);

% Extract domain parameters - define computational domain boundaries
x_min = cfg.domain.x_min;  % Left boundary of domain
x_max = cfg.domain.x_max;  % Right boundary of domain
y_min = cfg.domain.y_min;  % Bottom boundary of domain
y_max = cfg.domain.y_max;  % Top boundary of domain

%% 2) Generate geometry
G = build_geometry(cfg);

% Extract mesh data from geometry structure
xy = G.xy;                    % Interior velocity nodes
xy_s = G.xy_s;               % Interior pressure nodes
xt = G.xt;                   % Triangle connectivity

% Extract velocity grid boundaries
boundary_in = G.boundary_in;     % Inlet boundary nodes (V-grid)
boundary_out = G.boundary_out;   % Outlet boundary nodes (V-grid)
boundary_y = G.boundary_y;       % Wall boundary nodes (V-grid)
boundary_obs = G.boundary_obs;   % Obstacle boundary nodes (V-grid)

% Extract pressure grid boundaries
boundary_in_s = G.boundary_in_s;   % Inlet boundary nodes (P-grid)
boundary_out_s = G.boundary_out_s; % Outlet boundary nodes (P-grid)
boundary_y_s = G.boundary_y_s;     % Wall boundary nodes (P-grid)
boundary_obs_s = G.boundary_obs_s; % Obstacle boundary nodes (P-grid)

% Extract geometry-specific parameters (used for distance calculations)
if isfield(G, 'radius')
  radius = G.radius;           % Cylinder radius
  obs_char_length = radius;    % Characteristic length for distance calculations
elseif isfield(G, 'a') && isfield(G, 'b')
  a = G.a;                     % Ellipse semi-major axis
  b = G.b;                     % Ellipse semi-minor axis
  obs_char_length = min(a, b); % Use minimum axis as characteristic length
elseif isfield(G, 'rect_width') && isfield(G, 'rect_height')
  rect_width = G.rect_width;   % Rectangle width
  rect_height = G.rect_height; % Rectangle height
  obs_char_length = min(rect_width, rect_height) / 2; % Use half of smaller dimension
elseif isfield(G, 'chord_length') && isfield(G, 'naca_digits')
  chord_length = G.chord_length; % Airfoil chord length
  naca_digits = G.naca_digits;   % NACA 4-digit parameters
  % Calculate maximum thickness as characteristic length
  thickness_percent = naca_digits(3) * 10 + naca_digits(4); % Thickness as percentage
  max_thickness = (thickness_percent / 100) * chord_length;  % Maximum thickness in units
  obs_char_length = max_thickness / 2; % Use half of maximum thickness
elseif isfield(G, 'obstacles') && isfield(G, 'num_obstacles')
  % Multi-obstacle geometry: use minimum characteristic length of all obstacles
  obs_char_length = inf;
  for i = 1:G.num_obstacles
    obs = G.obstacles(i);
    switch lower(obs.type)
      case 'cylinder'
        char_len = obs.params.radius;
      case 'ellipse'
        char_len = min(obs.params.a, obs.params.b);
      case 'rectangle'
        char_len = min(obs.params.width, obs.params.height) / 2;
      case 'airfoil'
        thickness_percent = obs.params.naca_digits(3) * 10 + obs.params.naca_digits(4);
        max_thickness = (thickness_percent / 100) * obs.params.chord_length;
        char_len = max_thickness / 2;
      otherwise
        char_len = 0.5; % Default fallback
    end
    obs_char_length = min(obs_char_length, char_len);
  end
else
  error('Geometry structure missing expected parameters');
end

% Combine boundary nodes in specific order for compatibility
boundary_s = [boundary_y_s; boundary_out_s; boundary_in_s; boundary_obs_s];
boundary = [boundary_y; boundary_out; boundary_in; boundary_obs];

% Initial mesh visualization (only in non-CI mode)
if doPlot && ~isCI && ~isTest
  fprintf('Creating mesh visualization...\n');
  fig = figure('Name', 'Multi-Obstacle Mesh', 'Position', [100, 100, 1000, 600]);

  % Plot interior nodes (sample for performance if needed)
  if size(xy, 1) > 10000
    % Sample for very large meshes
    sample_rate = max(1, floor(size(xy, 1) / 8000));
    scatter(xy(1:sample_rate:end, 1), xy(1:sample_rate:end, 2), 4, 'k.', 'DisplayName', ...
            sprintf('Interior nodes (sampled %d/%d)', length(1:sample_rate:size(xy, 1)), size(xy, 1)));
  else
    % Show all interior nodes for reasonable mesh sizes
    scatter(xy(:, 1), xy(:, 2), 3, 'k.', 'DisplayName', sprintf('Interior nodes (%d)', size(xy, 1)));
  end
  hold on;

  % Plot pressure grid interior nodes (lighter color)
  if size(xy_s, 1) > 5000
    % Sample pressure nodes for very large meshes
    p_sample_rate = max(1, floor(size(xy_s, 1) / 4000));
    scatter(xy_s(1:p_sample_rate:end, 1), xy_s(1:p_sample_rate:end, 2), 2, [0.7 0.7 0.7], '.', ...
            'DisplayName', sprintf('Pressure nodes (sampled %d/%d)', length(1:p_sample_rate:size(xy_s, 1)), size(xy_s, 1)));
  else
    scatter(xy_s(:, 1), xy_s(:, 2), 2, [0.7 0.7 0.7], '.', ...
            'DisplayName', sprintf('Pressure nodes (%d)', size(xy_s, 1)));
  end

  % Plot boundary nodes with different colors and markers
  scatter(boundary_in(:, 1), boundary_in(:, 2), 25, 'b', 's', 'filled', 'DisplayName', sprintf('Inlet (%d)', length(boundary_in)));
  scatter(boundary_y(:, 1), boundary_y(:, 2), 25, 'r', '^', 'filled', 'DisplayName', sprintf('Walls (%d)', length(boundary_y)));
  scatter(boundary_out(:, 1), boundary_out(:, 2), 25, 'g', 'v', 'filled', 'DisplayName', sprintf('Outlet (%d)', length(boundary_out)));

  % Color-code obstacle boundaries if multi-obstacle
  if isfield(G, 'boundary_obs_ids') && isfield(G, 'num_obstacles')
    colors = {[1 0 1], [0 1 1], [1 1 0], [0.8 0.4 0], [0.6 0.2 0.8], [0.2 0.8 0.2], [0.8 0.2 0.2]};
    for obs_id = 1:G.num_obstacles
      obs_mask = G.boundary_obs_ids == obs_id;
      if any(obs_mask)
        color_idx = mod(obs_id - 1, length(colors)) + 1;
        color = colors{color_idx};
        scatter(G.boundary_obs_s(obs_mask, 1), G.boundary_obs_s(obs_mask, 2), ...
                35, color, 'o', 'filled', 'LineWidth', 1, 'MarkerEdgeColor', 'k', ...
                'DisplayName', sprintf('Obstacle %d (%d pts)', obs_id, sum(obs_mask)));
      end
    end
  else
    scatter(boundary_obs(:, 1), boundary_obs(:, 2), 35, [1 0 1], 'o', 'filled', ...
            'LineWidth', 1, 'MarkerEdgeColor', 'k', ...
            'DisplayName', sprintf('Obstacle (%d pts)', length(boundary_obs)));
  end

  axis equal;
  grid on;

  % Calculate total node counts
  total_v_nodes = size(xy, 1) + length(boundary_in) + length(boundary_y) + length(boundary_out) + length(boundary_obs);
  total_p_nodes = size(xy_s, 1) + length(G.boundary_y_s) + length(G.boundary_out_s) + length(G.boundary_in_s) + length(G.boundary_obs_s);
  if isfield(G, 'num_obstacles')
    num_obstacles = G.num_obstacles;
  else
    num_obstacles = 1; % Single obstacle case
  end

  title(sprintf('Multi-Obstacle Mesh: %d obstacles, V-grid: %d nodes, P-grid: %d nodes', ...
                num_obstacles, total_v_nodes, total_p_nodes));
  xlabel('x');
  ylabel('y');
  legend('Location', 'best', 'FontSize', 8);

  % Save figure and display
  try
    saveas(fig, 'multi_obstacle_mesh.png');
    fprintf('Mesh visualization saved as multi_obstacle_mesh.png\n');
  catch
    % Ignore save errors
  end

  % Force figure to front and pause briefly to ensure visibility
  figure(fig);
  drawnow;
  fprintf('Mesh visualization displayed. Continuing with simulation...\n');
  pause(1); % Brief pause to let user see the mesh
elseif isCI || isTest
  fprintf('Skipping mesh visualization (CI/test mode)\n');
end

% Combine interior and boundary nodes to create complete grids
xy1 = [xy; boundary];      % Complete velocity grid (V-grid): interior + boundary
xy1_s = [xy_s; boundary_s]; % Complete pressure grid (P-grid): interior + boundary

% Extract coordinates for convenience
x1 = xy1(:, 1);   % x-coordinates of velocity nodes
y1 = xy1(:, 2);   % y-coordinates of velocity nodes
x1_s = xy1_s(:, 1); % x-coordinates of pressure nodes
y1_s = xy1_s(:, 2); % y-coordinates of pressure nodes
x0 = xy(:, 1);    % x-coordinates of interior velocity nodes only
y0 = xy(:, 2);    % y-coordinates of interior velocity nodes only
x0_s = xy_s(:, 1); % x-coordinates of interior pressure nodes only
y0_s = xy_s(:, 2); % y-coordinates of interior pressure nodes only

% Define distance thresholds for special RBF treatment near boundaries
x_min_dist = cfg.distances.x_min;  % Distance threshold from left boundary
x_max_dist = cfg.distances.x_max;  % Distance threshold from right boundary
y_min_dist = cfg.distances.y_min;  % Distance threshold from bottom boundary
y_max_dist = cfg.distances.y_max;  % Distance threshold from top boundary
r_dist = cfg.mesh.refine_a2 * cfg.mesh.edge_multiplier;  % Distance threshold from cylinder

%% 3) Build stencils for RBF-FD method
S = build_stencils(xy1, xy1_s, xy, xy_s, G, cfg);

%% 4) Build pressure system and boundary conditions
P = build_pressure_system(G, xy_s, xy1_s, cfg);
L_inv_s = P.L_inv_s;
D0_21_x_obs = P.D0_21_x_obs;
D0_21_y_obs = P.D0_21_y_obs;

%% 5) Build inter-grid operators (P-grid ↔ V-grid)
[D0_21_x, D0_21_y, D0_12_x, D0_12_y] = build_intergrid_ops(G, xy, xy1, xy_s, xy1_s, S, cfg);

%% 6) Build velocity operators and boundary conditions
Vops = build_velocity_operators(G, xy, xy1, boundary_y, boundary_out, S, cfg, cfg.distances);

Dx = Vops.Dx;
Dy = Vops.Dy;
L0 = Vops.L0;
Dy_b_0 = Vops.Dy_b_0;
Dx_b_0 = Vops.Dx_b_0;
Dy_b = Vops.Dy_b;
Dy_b_1 = Vops.Dy_b_1;

%% 7) Precompute velocity system solvers
[L_u_inv, L_v_inv] = precompute_velocity_solvers(L0, Dy_b_0, Dx_b_0, xy, boundary_y, boundary_out, cfg);

%% 8) Setup time-stepping parameters
nu = cfg.simulation.viscosity;  % Kinematic viscosity (1/Reynolds number)
dt = cfg.simulation.time_step;  % Time step size

%% 8.1) Precompute mesh spacing stats for CFL-like diagnostics
h_min = NaN;
h_med = NaN;
if cfg.logging.enable
  try
    [~, D2] = knnsearch(xy, xy, 2);  % 2nd neighbor distance on interior V-grid
    if size(D2, 2) >= 2
      d2 = D2(:, 2);
      h_min = min(d2);
      h_med = median(d2);
      fprintf('[diag] spacing: h_min=%.3g, h_med=%.3g (V-grid interior)\n', h_min, h_med);
    end
  catch ME
    warning('spacing precomp failed: %s', ME.message);
  end
end

%% 9) Initialize simulation state
[W, p0] = init_state(xy1, xy1_s, boundary_obs, Nt);

% Define useful index lengths for boundary handling
L_B = length(boundary_obs) + length(boundary_in);   % Total special boundaries
L_B_y = length(boundary_y);   % Wall boundaries
L_W = length(xy);             % Interior velocity nodes
L_B_obs = length(boundary_obs);   % Obstacle boundary nodes

% Define pressure boundary indices for easy access
Idx_y_s = length(xy_s) + 1:length(xy_s) + length(boundary_y_s);    % Wall pressure indices
Idx_in_s = length(xy1_s) + 1 - length(boundary_obs_s) - length(boundary_in_s): ...
           length(xy1_s) - length(boundary_obs_s);  % Inlet pressure indices
Idx_out_s = length(xy_s) + 1 + length(boundary_y_s): ...
            length(xy_s) + length(boundary_out_s) + length(boundary_y_s);  % Outlet pressure indices

%% 10) Main time-stepping loop: Fractional Step Method
for j = 1:Nt
  if cfg.simulation.show_progress && ~isCI && ~isTest
    disp(['Time step j = ' num2str(j)]);
  end

  try
    % Use different schemes for startup vs main integration
    if j < 3
      % First few steps: use simple first-order scheme for stability
      W(:, j + 1) = W(:, j);  % Copy previous solution for startup
    else
      % Main fractional step algorithm implemented in NS_2d_fractional_step_PHS
      % This performs: 1) Advection-diffusion step with Adams-Bashforth + Crank-Nicolson
      %                2) Pressure correction step to enforce incompressibility
      %                3) Velocity correction using pressure gradient
      [W(:, j + 1), p0] = NS_2d_fractional_step_PHS(dt, nu, W(:, j - 1), W(:, j), ...
                                                    Dy, Dx, L_inv_s, L_u_inv, L_v_inv, ...
                                                    L0, L_B, L_B_obs, L_W, L_B_y, ...
                                                    length(boundary_s), D0_12_x, D0_12_y, ...
                                                    D0_21_x, D0_21_y, Dy_b, Dy_b_1, ...
                                                    D0_21_x_obs, D0_21_y_obs, p0, W(:, 1));
    end
  catch ME
    warning('Step %d threw error: %s', j, ME.message);
    if cfg.logging.snapshot_on_nan
      save_debug_snapshot(cfg, j, xy1, W, p0, dt, Dx, Dy, D0_12_x, D0_12_y, G);
    end
    rethrow(ME);
  end

  % Diagnostics every N steps
  if cfg.logging.enable && (mod(j, cfg.logging.step_frequency) == 0 || j <= 3)
    log_step_stats(j, W(:, j), W(:, j + 1), dt, Dx, Dy, D0_12_x, D0_12_y, xy1, h_min, h_med, cfg);
  end

  % Instability detection and snapshot
  if any(~isfinite(W(:, j + 1)))
    bad_idx = find(~isfinite(W(:, j + 1)), 1, 'first');
    [loc_str, pt] = classify_node_idx(bad_idx, L_W, boundary_y, boundary_out, boundary_in, boundary_obs, xy1);

    fprintf('\n*** SIMULATION INSTABILITY DETECTED ***\n');
    fprintf('Step: %d, Index: %d, Location: %s at (%.3g, %.3g)\n', j, bad_idx, loc_str, pt(1), pt(2));

    % Perform detailed instability analysis
    analyze_instability_cause(j, W, p0, bad_idx, loc_str, pt, dt, Dx, Dy, D0_12_x, D0_12_y, xy1, h_min, h_med, cfg);

    if cfg.logging.snapshot_on_nan
      save_debug_snapshot(cfg, j, xy1, W, p0, dt, Dx, Dy, D0_12_x, D0_12_y, G);
    end
    break
  end
end

% Extract final velocity field for potential continuation
W0 = W(:, end);

%% 11) Visualization of final results
visualize_final(cfg, doPlot, xy1, W, Nt, x_min, x_max, y_min, y_max, Dx, Dy);
