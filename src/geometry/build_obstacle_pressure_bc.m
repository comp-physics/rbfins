function [Dn1_b_s, Nearest_Idx_b_obs] = build_obstacle_pressure_bc(G, xy_s, xy1_s, cfg)
%BUILD_OBSTACLE_PRESSURE_BC Build pressure boundary condition operators for obstacle
%
% This function builds the normal derivative operator for Neumann pressure boundary
% conditions on the obstacle surface. For incompressible flow, the pressure normal
% derivative on solid boundaries satisfies: dp/dn = 0 (no acceleration into wall).
%
% INPUTS:
%   G       - Geometry structure from make_cylinder_geometry (or other geometry helper)
%   xy_s    - Interior pressure grid nodes
%   xy1_s   - Complete pressure grid nodes [xy_s; boundary_s]
%   cfg     - Configuration structure
%
% OUTPUTS:
%   Dn1_b_s           - Normal derivative operator for obstacle boundary condition
%   Nearest_Idx_b_obs - Nearest neighbor indices for obstacle boundary nodes

% Find nearest neighbors for obstacle boundary pressure nodes
[Nearest_Idx_b_obs] = nearest_interp(G.boundary_obs_s, xy_s, cfg.rbf.stencil_size_boundary_obstacle);

% Add boundary node indices to neighbor list
% The obstacle boundary nodes are at the end of the complete pressure grid
obs_start_idx = length(xy1_s) - length(G.boundary_obs_s) + 1;
obs_end_idx = length(xy1_s);
obs_indices = (obs_start_idx:obs_end_idx)';

Nearest_Idx_b_obs = [obs_indices, Nearest_Idx_b_obs];

% Generate RBF-FD weights for gradient operators at obstacle boundary
D = RBF_PHS_FD_all(G.boundary_obs_s, xy1_s, Nearest_Idx_b_obs, ...
    cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);

% Compute normal derivative: dp/dn = (nx*dp/dx + ny*dp/dy)
% Uses precomputed unit normal vectors from geometry helper
Dn1_b_s = G.obs_normals_s(:,1) .* D{1} + G.obs_normals_s(:,2) .* D{2};

end
