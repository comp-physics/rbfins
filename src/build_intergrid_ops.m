function [D0_21_x, D0_21_y, D0_12_x, D0_12_y] = build_intergrid_ops(G, xy, xy1, xy_s, xy1_s, S, cfg)
%BUILD_INTERGRID_OPS Generate interpolation operators between grids
%
% This function builds RBF-FD operators for interpolation and differentiation
% between the velocity grid (V-grid) and pressure grid (P-grid).
%
% INPUTS:
%   G     - Geometry structure
%   xy    - Interior velocity nodes
%   xy1   - Complete velocity grid (interior + boundary)
%   xy_s  - Interior pressure nodes
%   xy1_s - Complete pressure grid (interior + boundary)
%   S     - Stencil structure from build_stencils
%   cfg   - Configuration structure
%
% OUTPUTS:
%   D0_21_x - x-derivative operator (P-grid to V-grid)
%   D0_21_y - y-derivative operator (P-grid to V-grid)
%   D0_12_x - x-derivative operator (V-grid to P-grid)
%   D0_12_y - y-derivative operator (V-grid to P-grid)

% Differentiation operators: P-grid to V-grid
% Low-order PHS-RBFs and polynomials near boundaries
[D0_21_all] = RBF_PHS_FD_all(xy1(1:length(xy)+length(G.boundary_y)+length(G.boundary_out), :), ...
    xy1_s, S.Nearest_Idx_interp_21, cfg.rbf.order_interpolation_low, cfg.rbf.poly_degree_main, ...
    cfg.rbf.poly_degree_interpolation_low);

D0_21_x = D0_21_all{1};
D0_21_y = D0_21_all{2};

% Use precomputed indices for nodes far from all boundaries (high-order region)
Nearest_Idx_nc = find(G.idx_far_boundaries_V);
xy_nc = xy1(Nearest_Idx_nc, :);

[D0_21_all_nc] = RBF_PHS_FD_all(xy_nc, xy1_s, S.Nearest_Idx_interp_21(Nearest_Idx_nc, :), ...
    cfg.rbf.order_interpolation_high, cfg.rbf.order_interpolation_high_poly, ...
    cfg.rbf.poly_degree_interpolation_high);

D0_21_x(Nearest_Idx_nc, :) = D0_21_all_nc{1};
D0_21_y(Nearest_Idx_nc, :) = D0_21_all_nc{2};

% Differentiation operators: V-grid to P-grid
% Low-order PHS-RBFs and polynomials near boundaries
D0_12_all = RBF_PHS_FD_all(xy_s, xy1, S.Nearest_Idx_interp, cfg.rbf.order_boundary, ...
    cfg.rbf.order_interpolation_high_poly, cfg.rbf.derivative_order);

D0_12_x = D0_12_all{1};
D0_12_y = D0_12_all{2};

% Use precomputed indices for pressure nodes far from all boundaries (high-order region)
Nearest_Idx_nc = find(G.idx_far_boundaries_P);
xy_nc = xy1_s(Nearest_Idx_nc, :);

D0_12_all_nc = RBF_PHS_FD_all(xy_nc, xy1, S.Nearest_Idx_interp(Nearest_Idx_nc, :), ...
    cfg.rbf.order_interpolation_high, cfg.rbf.order_interpolation_high_poly, ...
    cfg.rbf.poly_degree_interpolation_high);

D0_12_x(Nearest_Idx_nc, :) = D0_12_all_nc{1};
D0_12_y(Nearest_Idx_nc, :) = D0_12_all_nc{2};

end
