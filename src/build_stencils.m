function S = build_stencils(xy1, xy1_s, xy, xy_s, G, cfg)
%BUILD_STENCILS Determine local stencils for RBF-FD method
%
% This function computes k-nearest neighbor stencils for both velocity and
% pressure grids, which are used for RBF-FD differentiation operators.
%
% INPUTS:
%   xy1   - Complete velocity grid (interior + boundary nodes)
%   xy1_s - Complete pressure grid (interior + boundary nodes)
%   xy    - Interior velocity nodes only
%   xy_s  - Interior pressure nodes only
%   G     - Geometry structure
%   cfg   - Configuration structure
%
% OUTPUTS:
%   S - Structure containing stencil information:
%       .Nearest_Idx         - k-nearest neighbors for V-grid
%       .Nearest_Idx_s       - k-nearest neighbors for P-grid
%       .Nearest_Idx_interp  - Interpolation stencils (P-grid to V-grid)
%       .Nearest_Idx_interp_21 - Interpolation stencils (V-grid to P-grid)

% For velocity nodes: find k nearest neighbors for each node
k = cfg.rbf.stencil_size_main; % Number of neighbors for main stencils
S.Nearest_Idx = nearest_interp(xy1, xy1, k); % k-nearest neighbors for V-grid

% For pressure nodes: find k nearest neighbors for each node
k = cfg.rbf.stencil_size_main;
S.Nearest_Idx_s = nearest_interp(xy1_s, xy1_s, k); % k-nearest neighbors for P-grid

% Interpolation stencils: P-grid to V-grid
S.Nearest_Idx_interp_21 = nearest_interp(xy1(1:length(xy)+length(G.boundary_y)+length(G.boundary_out), :), ...
xy1_s, cfg.rbf.stencil_size_boundary_outlet);

% Interpolation stencils: V-grid to P-grid
S.Nearest_Idx_interp = nearest_interp(xy_s, xy1, cfg.rbf.stencil_size_boundary_outlet);

end
