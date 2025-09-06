function Vops = build_velocity_operators(G, xy, xy1, boundary_y, boundary_out, S, cfg, distances)
%BUILD_VELOCITY_OPERATORS Generate velocity operators and boundary conditions
%
% This function builds RBF-FD operators for velocity field computation including:
% - Spatial derivatives (Dx, Dy) and Laplacian (L0) on velocity grid
% - Boundary condition operators for walls and outlets
% - Special treatment near obstacles and domain boundaries
%
% INPUTS:
%   G           - Geometry structure
%   xy          - Interior velocity nodes
%   xy1         - Complete velocity grid (interior + boundary)
%   boundary_y  - Wall boundary nodes
%   boundary_out - Outlet boundary nodes
%   S           - Stencil structure from build_stencils
%   cfg         - Configuration structure
%   distances   - Distance thresholds structure
%
% OUTPUTS:
%   Vops - Structure containing velocity operators:
%          .Dx, .Dy, .L0     - Spatial derivative and Laplacian operators
%          .Dy_b_0, .Dx_b_0  - Boundary condition operators
%          .Dy_b, .Dy_b_1    - Modified wall boundary operators
%          .Dx_b, .Dx_b_1    - Modified outlet boundary operators

% Extract distance thresholds
x_min_dist = distances.x_min;
x_max_dist = distances.x_max;
y_min_dist = distances.y_min;
y_max_dist = distances.y_max;

% Extract domain boundaries from config (assuming they're available)
x_min = cfg.domain.x_min;
x_max = cfg.domain.x_max;
y_min = cfg.domain.y_min;
y_max = cfg.domain.y_max;

% Boundary condition operators for V-grid
% Wall boundary conditions (du/dy = 0, dv/dy = 0)
[Nearest_Idx_b_y] = nearest_interp(boundary_y, xy, cfg.rbf.stencil_size_boundary_obstacle);
Nearest_Idx_b_y = [(length(xy)+1:length(xy)+length(boundary_y))' Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y, xy1, Nearest_Idx_b_y, cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Vops.Dy_b_0 = D{2};

Vops.Dy_b = Vops.Dy_b_0;
Vops.Dy_b_1 = diag(Vops.Dy_b(:,Nearest_Idx_b_y(:,1)));
Vops.Dy_b(:,Nearest_Idx_b_y(:,1)) = zeros(length(boundary_y),length(boundary_y));

% Outlet boundary conditions (du/dx = 0, dv/dx = 0)
[Nearest_Idx_b_out] = nearest_interp(boundary_out, xy, cfg.rbf.stencil_size_boundary_obstacle);
Nearest_Idx_b_out = [(length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out)))' Nearest_Idx_b_out];

D = RBF_PHS_FD_all(boundary_out, xy1, Nearest_Idx_b_out, cfg.rbf.order_boundary, cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Vops.Dx_b_0 = D{1};

Vops.Dx_b = Vops.Dx_b_0;
Vops.Dx_b_1 = diag(Vops.Dx_b(:,Nearest_Idx_b_out(:,1)));
Vops.Dx_b(:,Nearest_Idx_b_out(:,1)) = zeros(length(boundary_out),length(boundary_out));

% Generate main differentiation operators on velocity grid
D_all = RBF_PHS_FD_all(xy1, xy1, S.Nearest_Idx, cfg.rbf.order_main, cfg.rbf.poly_degree_main, cfg.rbf.laplacian_order);

Vops.Dx = D_all{1};
Vops.Dy = D_all{2};
Vops.L0 = D_all{3};

% Use precomputed indices for velocity nodes near obstacle (special RBF treatment)
Nearest_Idx_nb = find(G.idx_near_obs_V);
xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb, xy1, S.Nearest_Idx(Nearest_Idx_nb,:), cfg.rbf.order_near_obstacle, cfg.rbf.order_near_obstacle_poly, cfg.rbf.poly_degree_near_obstacle);

Vops.Dx(Nearest_Idx_nb,:) = D_all_nb{1};
Vops.Dy(Nearest_Idx_nb,:) = D_all_nb{2};
Vops.L0(Nearest_Idx_nb,:) = D_all_nb{3};

% Low-order PHS-RBFs and polynomials near domain boundaries
Nearest_Idx_nb = xy(:,1)>x_max-x_max_dist | xy(:,1)<x_min+x_min_dist*2 | xy(:,2)> y_max-y_max_dist | xy(:,2)<y_min+y_min_dist;
Idx = (1:length(xy))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);
xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb, xy1, S.Nearest_Idx(Nearest_Idx_nb,:), cfg.rbf.order_near_boundary, cfg.rbf.poly_degree_main, cfg.rbf.poly_degree_near_boundary);

Vops.Dx(Nearest_Idx_nb,:) = D_all_nb{1};
Vops.Dy(Nearest_Idx_nb,:) = D_all_nb{2};
Vops.L0(Nearest_Idx_nb,:) = D_all_nb{3};

end
