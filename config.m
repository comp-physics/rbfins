function config = config(geometry_type)
    % CONFIG - Configuration parameters for RBF Incompressible Navier-Stokes simulation
    %
    % This function returns a structure containing all configuration parameters
    % for flow simulations with different obstacle geometries.
    %
    % Usage:
    %   cfg = config();                    % Default geometry (cylinder)
    %   cfg = config('cylinder');          % Cylinder geometry
    %   cfg = config('ellipse');           % Ellipse geometry
    %   cfg = config('rectangle');         % Rectangle geometry
    %   cfg = config('airfoil');           % NACA airfoil geometry
    %   cfg = config('multi');             % Multiple obstacles geometry

    if nargin < 1
        geometry_type = 'multi';  % Default geometry
    end

    %% Domain Configuration
    config.domain.x_min = -8;
    config.domain.x_max = 24;
    config.domain.y_min = -8;
    config.domain.y_max = 8;

    %% Mesh Generation Parameters
    config.mesh.dist = 0.1;                   % Mesh spacing (adjusted per geometry below)
    config.mesh.boundary_eps = 0.002;         % Boundary tolerance
    config.mesh.refine_a1 = 0.05;             % Mesh refinement parameter A1 (near obstacle)
    config.mesh.refine_b1 = 0.08;             % Mesh refinement parameter B1 (near obstacle)
    config.mesh.refine_a2 = 0.05;             % Mesh refinement parameter A2 (wake region)
    config.mesh.refine_b2 = 0.08;             % Mesh refinement parameter B2 (wake region)
    config.mesh.edge_multiplier = 3;          % Multiplier for edge generation from triangles

    %% Simulation Parameters (set defaults first, then override per geometry)
    config.simulation.reynolds_number = 100;
    config.simulation.viscosity = 1 / config.simulation.reynolds_number;
    config.simulation.time_step = 1e-2;  % Default time step (may be overridden by geometry)
    config.simulation.num_time_steps = 2000;
    config.simulation.num_time_steps_ci = 20;
    config.simulation.random_seed = 42;  % Required for DistMesh reproducibility (uses rand() for rejection method)
    config.simulation.show_progress = true;  % Display time step progress (disabled in CI)

    %% Geometry Parameters - Set based on geometry type
    % Normalize geometry type to lowercase once
    geometry_type_lower = lower(geometry_type);
    config.geometry.type = geometry_type_lower;

    switch geometry_type_lower
        case 'cylinder'
            config.geometry.obstacle_radius = 0.5;    % Cylinder radius

        case 'ellipse'
            config.geometry.ellipse_a = 0.5;          % Ellipse semi-major axis (x-direction)
            config.geometry.ellipse_b = 0.4;          % Ellipse semi-minor axis (y-direction)

        case 'rectangle'
            config.geometry.rect_width = 1.0;         % Rectangle width (x-direction)
            config.geometry.rect_height = 0.8;        % Rectangle height (y-direction)
            config.geometry.rect_x_center = 0.0;      % Rectangle center X-coordinate
            config.geometry.rect_y_center = 0.0;      % Rectangle center Y-coordinate
            % Rectangle needs finer mesh for convergence
            config.mesh.dist = 0.05;

        case 'airfoil'
            % NACA 4-digit series parameters
            config.geometry.naca_digits = [0, 0, 1, 9];       % NACA airfoil
            config.geometry.chord_length = 2.3;               % Chord length
            config.geometry.angle_of_attack = -10;            % Angle of attack in degrees
            config.geometry.airfoil_x_center = -0.5;          % Airfoil center X-coordinate (leading edge, shifted left)
            config.geometry.airfoil_y_center = 0.0;           % Airfoil center Y-coordinate

            config.mesh.dist = 0.075;                          % Need a somewhat finer mesh distribution
            config.simulation.time_step = 5e-3;               % Airfoil needs smaller time step for stability

        case 'multi'
            % Multiple obstacles configuration
            % Default: two cylinders side-by-side vertically
            config.geometry.obstacles = [ ...
                                         struct('type', 'cylinder', 'center', [0, 2.5], 'params', struct('radius', 0.5)), ...
                                         struct('type', 'cylinder', 'center', [0, -2.5], 'params', struct('radius', 0.5)) ...
                                        ];

            config.mesh.dist = 0.07;                          % Need a somewhat finer mesh distribution
            config.simulation.time_step = 5e-3;

        otherwise
            error('Unknown geometry type: %s. Supported types: cylinder, ellipse, rectangle, airfoil, multi', geometry_type);
    end

    %% RBF-FD Algorithm Parameters
    % Main stencil sizes
    config.rbf.stencil_size_main = 35;
    config.rbf.stencil_size_boundary_obstacle = 20;
    config.rbf.stencil_size_boundary_wall = 15;
    config.rbf.stencil_size_boundary_outlet = 30;

    % RBF orders and polynomial degrees
    config.rbf.order_main = 28;
    config.rbf.poly_degree_main = 3;
    config.rbf.laplacian_order = 3;

    config.rbf.order_boundary = 8;
    config.rbf.poly_degree_boundary = 3;
    config.rbf.derivative_order = 1;

    config.rbf.order_interpolation_low = 12;
    config.rbf.poly_degree_interpolation_low = 2;

    config.rbf.order_interpolation_high = 25;
    config.rbf.order_interpolation_high_poly = 5;
    config.rbf.poly_degree_interpolation_high = 3;

    config.rbf.order_near_obstacle = 30;
    config.rbf.order_near_obstacle_poly = 7;
    config.rbf.poly_degree_near_obstacle = 4;

    config.rbf.order_near_boundary = 17;
    config.rbf.poly_degree_near_boundary = 2;

    %% Boundary Distance Thresholds for Special Treatment
    config.distances.x_min = 1;
    config.distances.x_max = 1;
    config.distances.y_min = 0.5;
    config.distances.y_max = 0.5;

    %% Numerical Scheme Coefficients
    config.schemes.adams_bashforth_current = 3 / 2;     % 3/2 coefficient for current time step
    config.schemes.adams_bashforth_previous = 1 / 2;    % 1/2 coefficient for previous time step
    config.schemes.crank_nicolson = 1 / 2;              % 1/2 coefficient for Crank-Nicolson

    %% Visualization Parameters
    config.visualization.scatter_size = 15;
    config.visualization.plot_tick_y = [-5, 0, 5];
    config.visualization.plot_tick_x = [-5, 0, 5, 10, 15];
    config.visualization.color_axis_range = 1e-0;

end
