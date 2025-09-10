function [D0_21_x_obs, D0_21_y_obs] = build_obstacle_grad_operators(G, xy1_s, cfg)
  %BUILD_OBSTACLE_GRAD_OPERATORS Build gradient interpolation operators for obstacle boundary
  %
  % This function builds RBF-FD gradient operators for interpolation from pressure
  % grid (P-grid) to velocity grid obstacle boundary nodes. These operators are used
  % for pressure-dependent boundary conditions on the obstacle surface.
  %
  % INPUTS:
  %   G       - Geometry structure from make_cylinder_geometry (or other geometry helper)
  %   xy1_s   - Complete pressure grid nodes [xy_s; boundary_s]
  %   cfg     - Configuration structure
  %
  % OUTPUTS:
  %   D0_21_x_obs - x-gradient operator for obstacle boundary (P-grid to V-grid)
  %   D0_21_y_obs - y-gradient operator for obstacle boundary (P-grid to V-grid)

  % Find nearest neighbors for interpolation from P-grid to obstacle boundary nodes
  [Nearest_Idx_interp_21_obs] = nearest_interp(G.boundary_obs, xy1_s, cfg.rbf.stencil_size_boundary_obstacle);

  % Build RBF-FD gradient operators for obstacle boundary
  [D0_21_all_obs] = RBF_PHS_FD_all(G.boundary_obs, xy1_s, Nearest_Idx_interp_21_obs, ...
                                   cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);

  % Extract x and y gradient operators
  D0_21_x_obs = D0_21_all_obs{1}; % x-gradient operator (dp/dx)
  D0_21_y_obs = D0_21_all_obs{2}; % y-gradient operator (dp/dy)

end
