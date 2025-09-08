function [L_u_inv, L_v_inv] = precompute_velocity_solvers(L0, Dy_b_0, Dx_b_0, xy, boundary_y, boundary_out, cfg)
    %PRECOMPUTE_VELOCITY_SOLVERS Create velocity system solvers with boundary conditions
    %
    % This function creates the implicit diffusion operators for the Crank-Nicolson
    % scheme and precomputes their LU decompositions for efficient solving during
    % time stepping.
    %
    % INPUTS:
    %   L0       - Laplacian operator on velocity grid
    %   Dy_b_0   - Wall boundary condition operator (du/dy = 0)
    %   Dx_b_0   - Outlet boundary condition operator (du/dx = 0)
    %   xy       - Interior velocity nodes
    %   boundary_y  - Wall boundary nodes
    %   boundary_out - Outlet boundary nodes
    %   cfg      - Configuration structure
    %
    % OUTPUTS:
    %   L_u_inv - Function handle for solving u-velocity system
    %   L_v_inv - Function handle for solving v-velocity system

    % Extract simulation parameters
    nu = cfg.simulation.viscosity; % Kinematic viscosity (1/Reynolds number)
    dt = cfg.simulation.time_step; % Time step size

    % Compute this from the grid size
    xy1_size = size(L0, 1); % Total grid size from Laplacian operator
    boundary_size = xy1_size - length(xy);

    % Setup Laplacian operator for velocity diffusion
    % Zero out Laplacian on boundary nodes (boundary conditions applied separately)
    L = L0;
    L(length(xy) + 1:end, :) = zeros(boundary_size, xy1_size);

    % Create implicit diffusion operator for Crank-Nicolson scheme
    % I - (dt*nu/2)*del^2 for implicit half of viscous terms
    L_I = speye(xy1_size) - dt * nu * L * cfg.schemes.crank_nicolson;

    % Create separate operators for u and v velocity components
    L_u = L_I; % Copy base diffusion operator for u-velocity
    L_v = L_I; % Copy base diffusion operator for v-velocity

    % Apply boundary conditions for u-velocity
    % Wall boundary: du/dy = 0
    wall_start = length(xy) + 1;
    wall_end = length(xy) + length(boundary_y);
    L_u(wall_start:wall_end, :) = Dy_b_0;

    % Outlet boundary: du/dx = 0
    outlet_start = length(xy) + length(boundary_y) + 1;
    outlet_end = length(xy) + length(boundary_y) + length(boundary_out);
    L_u(outlet_start:outlet_end, :) = Dx_b_0;

    % Apply boundary conditions for v-velocity
    % Outlet boundary: dv/dx = 0
    L_v(outlet_start:outlet_end, :) = Dx_b_0;

    % Note: Wall and inlet boundaries have Dirichlet conditions (v=0) applied directly
    % Obstacle boundary has no-slip condition (u=v=0) applied directly

    % Precompute LU decompositions for efficient velocity solves
    [LL, UU, pp, qq, rr] = lu(L_u);
    L_u_inv = @(v) (qq * (UU \ (LL \ (pp * (rr \ (v)))))); % u-velocity solver

    [LL, UU, pp, qq, rr] = lu(L_v);
    L_v_inv = @(v) (qq * (UU \ (LL \ (pp * (rr \ (v)))))); % v-velocity solver

end
