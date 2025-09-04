function config = config()
% CONFIG - Configuration parameters for RBF Incompressible Navier-Stokes simulation
%
% This function returns a structure containing all configuration parameters
% for the cylinder flow simulation. Modify values here to change simulation
% parameters without touching the main code.
%
% Usage: cfg = config();

%% Domain Configuration
config.domain.x_min = -8;
config.domain.x_max = 24;
config.domain.y_min = -8;
config.domain.y_max = 8;

%% Mesh Generation Parameters
config.mesh.dist = 0.1;                    % Mesh spacing
config.mesh.cylinder_radius = 0.5;         % Cylinder radius
config.mesh.boundary_eps = 0.002;          % Boundary tolerance
config.mesh.refine_a1 = 0.05;             % Mesh refinement parameter A1
config.mesh.refine_b1 = 0.08;             % Mesh refinement parameter B1
config.mesh.refine_a2 = 0.05;             % Mesh refinement parameter A2
config.mesh.refine_b2 = 0.08;             % Mesh refinement parameter B2
config.mesh.edge_multiplier = 3;          % Multiplier for edge generation from triangles

%% RBF-FD Algorithm Parameters
% Main stencil sizes
config.rbf.stencil_size_main = 35;
config.rbf.stencil_size_boundary_cylinder = 20;
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

config.rbf.order_near_cylinder = 30;
config.rbf.order_near_cylinder_poly = 7;
config.rbf.poly_degree_near_cylinder = 4;

config.rbf.order_near_boundary = 17;
config.rbf.poly_degree_near_boundary = 2;

%% Simulation Parameters
config.simulation.reynolds_number = 100;
config.simulation.viscosity = 1/config.simulation.reynolds_number;
config.simulation.time_step = 1e-2;
config.simulation.num_time_steps = 5000;
config.simulation.num_time_steps_ci = 20;
config.simulation.random_seed = 42;  % Required for DistMesh reproducibility (uses rand() for rejection method)

%% Distance Thresholds for Special Treatment
config.distances.x_min = 1;
config.distances.x_max = 1;
config.distances.y_min = 0.5;
config.distances.y_max = 0.5;

%% Numerical Scheme Coefficients
config.schemes.adams_bashforth_current = 3/2;     % 3/2 coefficient for current time step
config.schemes.adams_bashforth_previous = 1/2;    % 1/2 coefficient for previous time step
config.schemes.crank_nicolson = 1/2;              % 1/2 coefficient for Crank-Nicolson

%% Visualization Parameters
config.visualization.scatter_size = 15;
config.visualization.plot_tick_y = [-5, 0, 5];
config.visualization.plot_tick_x = [-5, 0, 5, 10, 15];
config.visualization.color_axis_range = 1e-0;

end
