function P = build_pressure_system(G, xy_s, xy1_s, cfg)
%BUILD_PRESSURE_SYSTEM Generate pressure system and boundary conditions
%
% This function builds the complete pressure Poisson system including:
% - Laplacian operator for interior nodes
% - Boundary condition operators for all boundaries
% - Assembled system matrix with regularization
% - Precomputed LU decomposition for efficient solving
%
% INPUTS:
%   G     - Geometry structure
%   xy_s  - Interior pressure nodes
%   xy1_s - Complete pressure grid (interior + boundary)
%   cfg   - Configuration structure
%
% OUTPUTS:
%   P - Structure containing pressure system components:
%       .L_inv_s        - Pressure solver function handle
%       .D0_21_x_obs    - Obstacle x-gradient operator
%       .D0_21_y_obs    - Obstacle y-gradient operator
%       .L_s            - Laplacian operator

% Extract boundary data from geometry
boundary_y_s = G.boundary_y_s;
boundary_out_s = G.boundary_out_s;
boundary_in_s = G.boundary_in_s;
boundary_obs_s = G.boundary_obs_s;
boundary_s = [boundary_y_s; boundary_out_s; boundary_in_s; boundary_obs_s];

% Build k-nearest neighbor stencils for pressure grid
k = cfg.rbf.stencil_size_main;
Nearest_Idx_s = nearest_interp(xy1_s, xy1_s, k);

% Generate Laplacian operator for pressure Poisson equation (P-grid to P-grid)
[D_s_all] = RBF_PHS_FD_all(xy_s, xy1_s, Nearest_Idx_s(1:length(xy_s), :), cfg.rbf.order_main, cfg.rbf.poly_degree_main, cfg.rbf.laplacian_order);
P.L_s = D_s_all{3}; % Extract Laplacian operator (del^2 p = d^2p/dx^2 + d^2p/dy^2)

% Generate boundary condition operators for pressure system
% Obstacle boundary: Neumann BC (dp/dn = 0, normal derivative = 0)
[Dn1_b_s, ~] = build_obstacle_pressure_bc(G, xy_s, xy1_s, cfg);

% Wall boundaries (top/bottom): Neumann BC (dp/dy = 0)
[Nearest_Idx_b_y] = nearest_interp(boundary_y_s, [xy_s; xy1_s(length(xy_s)+length(boundary_y_s)+1, :)], cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_y = [(length(xy_s) + 1:length(xy_s) + length(boundary_y_s))', Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y_s, xy1_s, Nearest_Idx_b_y, cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dy_b_s = D{2}; % y-derivative operator for wall boundary condition

% Inlet boundary: Neumann BC (dp/dx = 0)
[Nearest_Idx_b_x] = nearest_interp(boundary_in_s, xy1_s(1:length(xy_s)+length(boundary_y_s), :), cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_x = [(length(xy1_s) + 1 - length(boundary_obs_s) - length(boundary_in_s):length(xy1_s) - length(boundary_obs_s))', Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_in_s, xy1_s, Nearest_Idx_b_x, cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_in_s = D{1}; % x-derivative operator for inlet boundary condition

% Outlet boundary: Neumann BC (dp/dx = 0)
[Nearest_Idx_b_x] = nearest_interp(boundary_out_s, xy1_s(1:length(xy_s)+length(boundary_y_s), :), cfg.rbf.stencil_size_boundary_outlet);
Nearest_Idx_b_x = [(length(xy_s) + 1 + length(boundary_y_s):length(xy_s) + length(boundary_y_s) + length(boundary_out_s))', Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_out_s, xy1_s, Nearest_Idx_b_x, cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_out_s = D{1}; % x-derivative operator for outlet boundary condition

% Assemble pressure boundary condition matrices
% Initialize boundary condition matrices (each row corresponds to one boundary node)
L_bc_c = zeros(length(boundary_s), length(xy1_s)); % Obstacle BC matrix
L_bc_out = zeros(length(boundary_s), length(xy1_s)); % Outlet BC matrix
L_bc_y = zeros(length(boundary_s), length(xy1_s)); % Wall BC matrix
L_bc_in = zeros(length(boundary_s), length(xy1_s)); % Inlet BC matrix

% Fill obstacle boundary condition matrix (dp/dn = 0)
Idx_boundary_obs = length(xy1_s) - length(boundary_obs_s) + (1:length(boundary_obs_s));
L_bc_c(Idx_boundary_obs-length(xy_s), :) = Dn1_b_s;

% Fill outlet boundary condition matrix (dp/dx = 0)
Idx_boundary_out = length(xy_s) + length(boundary_y_s) + (1:length(boundary_out_s));
L_bc_out(Idx_boundary_out-length(xy_s), :) = Dx_out_s;

% Fill wall boundary condition matrix (dp/dy = 0)
Idx_boundary_y = length(xy_s) + (1:length(boundary_y_s));
L_bc_y(Idx_boundary_y-length(xy_s), :) = Dy_b_s;

% Fill inlet boundary condition matrix (dp/dx = 0)
Idx_boundary_in = length(xy1_s) - length(boundary_obs_s) - length(boundary_in_s) + (1:length(boundary_in_s));
L_bc_in(Idx_boundary_in-length(xy_s), :) = Dx_in_s;

% Assemble complete pressure system: [Laplacian at interior; BCs at boundary]
L1 = [P.L_s(1:length(xy_s), :); L_bc_c + L_bc_out + L_bc_y + L_bc_in];

% Add regularization to fix pressure datum (pressure is defined up to constant)
% This adds the constraint sum(p) = 0 to make the system uniquely solvable
L1 = [L1, [ones(length(xy_s), 1); zeros(length(boundary_s), 1)]; [ones(1, length(xy_s)), zeros(1, length(boundary_s))], 0];

% Precompute LU decomposition for efficient pressure solves
[LL, UU, pp, qq, rr] = lu(L1);
P.L_inv_s = @(v) (qq * (UU \ (LL \ (pp * (rr \ (v)))))); % Pressure solver function

% Generate obstacle boundary gradient operators for pressure-dependent boundary conditions
[P.D0_21_x_obs, P.D0_21_y_obs] = build_obstacle_grad_operators(G, xy1_s, cfg);

end
