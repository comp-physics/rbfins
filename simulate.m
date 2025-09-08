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
else
    error('Geometry structure missing expected parameters');
end

% Combine boundary nodes in specific order for compatibility
boundary_s = [boundary_y_s; boundary_out_s; boundary_in_s; boundary_obs_s];
boundary = [boundary_y; boundary_out; boundary_in; boundary_obs];

% Initial mesh visualization
if doPlot
    figure;
    scatter(xy(:, 1), xy(:, 2), 'k.');
    hold on;
    axis square;
    scatter(boundary_in(:, 1), boundary_in(:, 2), 'b+');
    scatter(boundary_y(:, 1), boundary_y(:, 2), 'r+');
    scatter(boundary_out(:, 1), boundary_out(:, 2), 'b+');
    scatter(boundary_obs(:, 1), boundary_obs(:, 2), 'm+');
    axis equal;
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

    % Check for numerical instability
    if isnan(W(1, j + 1))
        warning('Simulation became unstable (NaN detected). Stopping at time step %d', j);
        break
    end
end

% Extract final velocity field for potential continuation
W0 = W(:, end);

%% 11) Visualization of final results
visualize_final(cfg, doPlot, xy1, W, Nt, x_min, x_max, y_min, y_max);
