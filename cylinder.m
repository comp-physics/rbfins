clc;
clear;

% Load Configuration
cfg = config();

% Check if running in CI environment or headless mode
isCI = strcmpi(getenv('CI'), 'true');
isTest = strcmpi(getenv('MATLAB_TEST'), 'true');

% Completely disable plotting for CI/test environments
doPlot = ~isCI && ~isTest;

if ~doPlot
    % Aggressively disable all graphics for tests
    set(0, 'DefaultFigureVisible', 'off');
    fprintf('Plotting completely disabled for testing/CI\n');
    hasDisplay = false;
else
    % Only check display if we might want to plot
    try
        % Try to create a figure to test if display is available
        f = figure('Visible', 'off');
        close(f);
        hasDisplay = true;
    catch
        hasDisplay = false;
        doPlot = false;
    end
end

% Set random seed for DistMesh reproducibility (DistMesh uses rand() for rejection method)
if isfield(cfg.simulation, 'random_seed') && ~isempty(cfg.simulation.random_seed)
    rng(cfg.simulation.random_seed);
    fprintf('Random seed set to %d for DistMesh reproducibility\n', cfg.simulation.random_seed);
end

% Ensure DistMesh is available
if ~exist('distmesh2d', 'file')
    % Try to add distmesh from the standard location
    if exist('distmesh', 'dir')
        addpath('distmesh');
    else
        error('distmesh directory not found. Please run setup_paths() or add distmesh to your path manually.');
    end
end

% Extract domain parameters - define computational domain boundaries
x_min = cfg.domain.x_min;  % Left boundary of domain
x_max = cfg.domain.x_max;  % Right boundary of domain  
y_min = cfg.domain.y_min;  % Bottom boundary of domain
y_max = cfg.domain.y_max;  % Top boundary of domain

% Extract mesh parameters for DistMesh generation
dist = cfg.mesh.dist;           % Target edge length for mesh
radius = cfg.mesh.cylinder_radius;  % Radius of circular cylinder
eps = cfg.mesh.boundary_eps;   % Tolerance for boundary detection
a1 = cfg.mesh.refine_a1;       % Refinement parameter 1 near cylinder
b1 = cfg.mesh.refine_b1;       % Refinement parameter 2 near cylinder
a2 = cfg.mesh.refine_a2;       % Refinement parameter 3 for wake region
b2 = cfg.mesh.refine_b2;       % Refinement parameter 4 for wake region

% Define distance functions for DistMesh
% fd: Signed distance function for domain geometry (rectangle minus circle)
fd = @(p) ddiff(drectangle(p,x_min, x_max,y_min,y_max),dcircle(p,0,0,radius));

% fd1: Edge length function for local mesh refinement
% Provides smaller elements near cylinder and in wake region behind it
fd1 = @(p) min(a1+b1*abs(dcircle(p,0,0,radius)),a2+b2*abs(dpoly(p,[radius 0; x_max 0])));

% Define fixed corner points to maintain rectangular domain shape
fix = [x_min,y_min; x_min,y_max; x_max, y_max; x_max, y_min];

% Generate triangular mesh using DistMesh algorithm
% xy_s: nodes for pressure (P-grid), xt: triangle connectivity
[xy_s,xt] = distmesh2d(fd,fd1,dist, [x_min,y_min; x_max, y_max],fix );

% Store mesh data
nodes = {xy_s,xt};
xy_s = nodes{1}; % Pressure grid nodes (P-grid)
xt = nodes{2};   % Triangle connectivity

% Generate velocity grid nodes (V-grid) by adding edge midpoints
% This creates a staggered grid arrangement for better pressure-velocity coupling
xy = zeros(length(xt)*cfg.mesh.edge_multiplier,2);
for j = 1:length(xt)
    % Add midpoint of each triangle edge to create velocity nodes
    xy((j-1)*cfg.mesh.edge_multiplier+1,:) = (xy_s(xt(j,1),:)+xy_s(xt(j,2),:))/2;  % Edge 1-2
    xy((j-1)*cfg.mesh.edge_multiplier+2,:) = (xy_s(xt(j,1),:)+xy_s(xt(j,3),:))/2;  % Edge 1-3
    xy((j-1)*cfg.mesh.edge_multiplier+3,:) = (xy_s(xt(j,2),:)+xy_s(xt(j,3),:))/2;  % Edge 2-3
end
xy = unique(xy,'rows');  % Remove duplicate nodes

% Classify and extract boundary nodes for pressure grid (P-grid)
% Remove corner nodes first to avoid duplicate boundary classification
idx_corners = (xy_s(:,1)<x_min+eps & xy_s(:,2)<y_min+eps) | ...
              (xy_s(:,1)<x_min+eps & xy_s(:,2)>y_max-eps) | ...
              (xy_s(:,1)>x_max-eps & xy_s(:,2)>y_max-eps) | ...
              (xy_s(:,1)>x_max-eps & xy_s(:,2)<y_min+eps);
xy_s(idx_corners,:) = [];

% Extract inlet boundary nodes (left wall, x = x_min)
idx_b_in = xy_s(:,1)<x_min+eps;
boundary_in_s = xy_s(idx_b_in,:);
xy_s(idx_b_in,:) = [];

% Extract wall boundary nodes (top and bottom walls, y = y_min/y_max)
idx_b_y = xy_s(:,2)>y_max-eps | xy_s(:,2)<y_min+eps;
boundary_y_s = xy_s(idx_b_y,:);
xy_s(idx_b_y,:) = [];

% Extract outlet boundary nodes (right wall, x = x_max)
idx_b_out = xy_s(:,1)>x_max-eps;
boundary_out_s = xy_s(idx_b_out,:);
xy_s(idx_b_out,:) = [];

% Extract cylinder boundary nodes (circular boundary)
idx_b_c = xy_s(:,1).^2+xy_s(:,2).^2<(radius+eps)^2;
boundary_c_s = xy_s(idx_b_c,:);
xy_s(idx_b_c,:) = [];

% Combine all pressure boundary nodes in specific order
boundary_s = [boundary_y_s; boundary_out_s; boundary_in_s; boundary_c_s];

% Classify and extract boundary nodes for velocity grid (V-grid)
% Extract inlet boundary nodes
idx_b_in = xy(:,1)<x_min+eps;
boundary_in = xy(idx_b_in,:);
xy(idx_b_in,:) = [];

% Extract wall boundary nodes (top and bottom)
idx_b_y = xy(:,2)>y_max-eps | xy(:,2)<y_min+eps;
boundary_y = xy(idx_b_y,:);
xy(idx_b_y,:) = [];

% Extract outlet boundary nodes
idx_b_out = xy(:,1)>x_max-eps;
boundary_out = xy(idx_b_out,:);
xy(idx_b_out,:) = [];

% Extract cylinder boundary nodes
idx_b_c = xy(:,1).^2+xy(:,2).^2<(radius+eps)^2;
boundary_c = xy(idx_b_c,:);
xy(idx_b_c,:) = [];

% Combine all velocity boundary nodes in same order as pressure grid
boundary = [boundary_y; boundary_out; boundary_in; boundary_c];

if doPlot
    figure;

    scatter(xy(:,1),xy(:,2),'k.'); hold on; axis square;
    scatter(boundary_in(:,1),boundary_in(:,2),'b+');
    scatter(boundary_y(:,1),boundary_y(:,2),'r+');
    scatter(boundary_out(:,1),boundary_out(:,2),'b+');
    scatter(boundary_c(:,1),boundary_c(:,2),'m+');

    axis equal;
end

% Combine interior and boundary nodes to create complete grids
xy1 = [xy; boundary];      % Complete velocity grid (V-grid): interior + boundary
xy1_s = [xy_s; boundary_s]; % Complete pressure grid (P-grid): interior + boundary

% Extract coordinates for convenience
x1 = xy1(:,1);   % x-coordinates of velocity nodes
y1 = xy1(:,2);   % y-coordinates of velocity nodes
x1_s = xy1_s(:,1); % x-coordinates of pressure nodes
y1_s = xy1_s(:,2); % y-coordinates of pressure nodes
x0 = xy(:,1);    % x-coordinates of interior velocity nodes only
y0 = xy(:,2);    % y-coordinates of interior velocity nodes only
x0_s = xy_s(:,1); % x-coordinates of interior pressure nodes only
y0_s = xy_s(:,2); % y-coordinates of interior pressure nodes only

% Define distance thresholds for special RBF treatment near boundaries
x_min_dist = cfg.distances.x_min;  % Distance threshold from left boundary
x_max_dist = cfg.distances.x_max;  % Distance threshold from right boundary
y_min_dist = cfg.distances.y_min;  % Distance threshold from bottom boundary
y_max_dist = cfg.distances.y_max;  % Distance threshold from top boundary
r_dist = a2*cfg.mesh.edge_multiplier;  % Distance threshold from cylinder

%% Determine local stencils for RBF-FD method
% For velocity nodes: find k nearest neighbors for each node
k = cfg.rbf.stencil_size_main;  % Number of neighbors for main stencils
[Nearest_Idx] = nearest_interp(xy1,xy1,k);  % k-nearest neighbors for V-grid

% For pressure nodes: find k nearest neighbors for each node
k = cfg.rbf.stencil_size_main;
[Nearest_Idx_s] = nearest_interp(xy1_s,xy1_s,k);  % k-nearest neighbors for P-grid

%% Generate Laplacian operator for pressure Poisson equation (P-grid to P-grid)
% Generate RBF-FD weights for Laplacian operator at interior pressure nodes
[D_s_all] = RBF_PHS_FD_all(xy_s,xy1_s,Nearest_Idx_s(1:length(xy_s),:),cfg.rbf.order_main,cfg.rbf.poly_degree_main,cfg.rbf.laplacian_order);
L_s = D_s_all{3};  % Extract Laplacian operator (del^2 p = d^2p/dx^2 + d^2p/dy^2)

%% Generate boundary condition operators for pressure system
% Cylinder boundary: Neumann BC (dp/dn = 0, normal derivative = 0)
[Nearest_Idx_b_c] = nearest_interp(boundary_c_s,xy_s,cfg.rbf.stencil_size_boundary_cylinder);
Nearest_Idx_b_c = [(length(xy1_s)-length(boundary_c_s)+1:length(xy1_s))' Nearest_Idx_b_c];

% Generate RBF-FD weights for gradient operators at cylinder boundary
D = RBF_PHS_FD_all(boundary_c_s,xy1_s,Nearest_Idx_b_c, cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
% Compute normal derivative: dp/dn = (nx*dp/dx + ny*dp/dy) where n = [x,y]/r is unit normal
Dn1_b_s = ((boundary_c_s(:,1)).*D{1}+(boundary_c_s(:,2)).*D{2})./radius;

% Wall boundaries (top/bottom): Neumann BC (dp/dy = 0)
[Nearest_Idx_b_y] = nearest_interp(boundary_y_s,[xy_s; xy1_s(length(xy_s)+length(boundary_y_s)+1,:)],cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_y = [(length(xy_s)+1:length(xy_s)+length(boundary_y_s))' Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y_s,xy1_s,Nearest_Idx_b_y,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dy_b_s = D{2};  % y-derivative operator for wall boundary condition

% Inlet boundary: Neumann BC (dp/dx = 0)
[Nearest_Idx_b_x] = nearest_interp(boundary_in_s,xy1_s(1:length(xy_s)+length(boundary_y_s),:),cfg.rbf.stencil_size_boundary_wall);
Nearest_Idx_b_x = [(length(xy1_s)+1-length(boundary_c_s)-length(boundary_in_s):length(xy1_s)-length(boundary_c_s))' Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_in_s,xy1_s,Nearest_Idx_b_x,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dx_in_s = D{1};  % x-derivative operator for inlet boundary condition

% Outlet boundary: Neumann BC (dp/dx = 0)
[Nearest_Idx_b_x] = nearest_interp(boundary_out_s,xy1_s(1:length(xy_s)+length(boundary_y_s),:),cfg.rbf.stencil_size_boundary_outlet);
Nearest_Idx_b_x = [(length(xy_s)+1+length(boundary_y_s):length(xy_s)+length(boundary_y_s)+length(boundary_out_s))' Nearest_Idx_b_x];

D = RBF_PHS_FD_all(boundary_out_s,xy1_s,Nearest_Idx_b_x, cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_out_s = D{1};  % x-derivative operator for outlet boundary condition

%% Assemble pressure boundary condition matrices
% Initialize boundary condition matrices (each row corresponds to one boundary node)
L_bc_c = zeros(length(boundary_s), length(xy1_s));   % Cylinder BC matrix
L_bc_out = zeros(length(boundary_s), length(xy1_s)); % Outlet BC matrix  
L_bc_y = zeros(length(boundary_s), length(xy1_s));   % Wall BC matrix
L_bc_in = zeros(length(boundary_s), length(xy1_s));  % Inlet BC matrix

% Fill cylinder boundary condition matrix (dp/dn = 0)
Idx_boundary_c = length(xy1_s)-length(boundary_c_s)+[1:length(boundary_c_s)];
L_bc_c(Idx_boundary_c-length(xy_s),:) = Dn1_b_s;

% Fill outlet boundary condition matrix (dp/dx = 0)
Idx_boundary_out = length(xy_s) + length(boundary_y_s) + [1:length(boundary_out_s)];
L_bc_out(Idx_boundary_out-length(xy_s),:) = Dx_out_s;

% Fill wall boundary condition matrix (dp/dy = 0)
Idx_boundary_y = length(xy_s) + [1:length(boundary_y_s)];
L_bc_y(Idx_boundary_y-length(xy_s),:) = Dy_b_s;

% Fill inlet boundary condition matrix (dp/dx = 0)
Idx_boundary_in = length(xy1_s) - length(boundary_c_s) - length(boundary_in_s) + [1:length(boundary_in_s)];
L_bc_in(Idx_boundary_in-length(xy_s),:) = Dx_in_s;

% Assemble complete pressure system: [Laplacian at interior; BCs at boundary]
L1 = [L_s(1:length(xy_s),:); L_bc_c+L_bc_out+L_bc_y+L_bc_in];

% Add regularization to fix pressure datum (pressure is defined up to constant)
% This adds the constraint sum(p) = 0 to make the system uniquely solvable
L1 = [L1 [ones(length(xy_s),1); zeros(length(boundary_s),1)]; [ones(1,length(xy_s)) zeros(1,length(boundary_s))] 0];

% Precompute LU decomposition for efficient pressure solves
[LL,UU,pp,qq,rr] = lu(L1);
L_inv_s = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));  % Pressure solver function

%% Generate interpolation operators between grids (P-grid to V-grid)
[Nearest_Idx_interp_21_c] = nearest_interp(boundary_c,xy1_s,cfg.rbf.stencil_size_boundary_outlet);
[D0_21_all_c] = RBF_PHS_FD_all(boundary_c,xy1_s,Nearest_Idx_interp_21_c,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);

D0_21_x_c =  D0_21_all_c{1};
D0_21_y_c =  D0_21_all_c{2};

%%  differentiation : P-grid to V-grid
[Nearest_Idx_interp_21] = nearest_interp(xy1(1:length(xy)+length(boundary_y)+length(boundary_out),:),xy1_s,cfg.rbf.stencil_size_boundary_outlet);

% low-order PHS-RBFs and polynomials near boundaries
[D0_21_all] = RBF_PHS_FD_all(xy1(1:length(xy)+length(boundary_y)+length(boundary_out),:),xy1_s,Nearest_Idx_interp_21,cfg.rbf.order_interpolation_low,cfg.rbf.poly_degree_main,cfg.rbf.poly_degree_interpolation_low);

D0_21_x = D0_21_all{1};
D0_21_y = D0_21_all{2};

Idx = (1:length(xy))';

Nearest_Idx_nc = xy(:,1).^2+xy(:,2).^2>(radius+r_dist)^2  & xy(:,1)>x_min+x_min_dist & xy(:,2)<y_max-y_max_dist & xy(:,2)>y_min+y_min_dist & xy(:,1)< x_max-x_max_dist  ;
Nearest_Idx_nc = Idx(Nearest_Idx_nc);
xy_nc = xy1(Nearest_Idx_nc,:);

[D0_21_all_nc] = RBF_PHS_FD_all(xy_nc,xy1_s,Nearest_Idx_interp_21(Nearest_Idx_nc,:),cfg.rbf.order_interpolation_high,cfg.rbf.order_interpolation_high_poly,cfg.rbf.poly_degree_interpolation_high);

D0_21_x(Nearest_Idx_nc,:) =  D0_21_all_nc{1};
D0_21_y(Nearest_Idx_nc,:) =  D0_21_all_nc{2};

% Differentiation operators: V-grid to P-grid
[Nearest_Idx_interp] = nearest_interp(xy_s,xy1,cfg.rbf.stencil_size_boundary_outlet);

% low-order PHS-RBFs and polynomials near boundaries
D0_12_all = RBF_PHS_FD_all(xy_s,xy1,Nearest_Idx_interp,cfg.rbf.order_boundary,cfg.rbf.order_interpolation_high_poly,cfg.rbf.derivative_order);

D0_12_x =  D0_12_all{1};
D0_12_y =  D0_12_all{2};

Idx = (1:length(xy_s))';

Nearest_Idx_nc = xy_s(:,1).^2+xy_s(:,2).^2>(radius+r_dist)^2  & xy_s(:,1)>x_min+x_min_dist & xy_s(:,2)<y_max-y_max_dist & xy_s(:,2)>y_min+y_min_dist & xy_s(:,1)< x_max-x_max_dist;

Nearest_Idx_nc = Idx(Nearest_Idx_nc);
xy_nc = xy1_s(Nearest_Idx_nc,:);

D0_12_all_nc = RBF_PHS_FD_all(xy_nc,xy1,Nearest_Idx_interp(Nearest_Idx_nc,:),cfg.rbf.order_interpolation_high,cfg.rbf.order_interpolation_high_poly,cfg.rbf.poly_degree_interpolation_high);

D0_12_x(Nearest_Idx_nc,:)  =  D0_12_all_nc{1};
D0_12_y(Nearest_Idx_nc,:)  =  D0_12_all_nc{2};

%   B.C.s for V-grid
[Nearest_Idx_b_y] = nearest_interp(boundary_y,xy,cfg.rbf.stencil_size_boundary_cylinder);
Nearest_Idx_b_y = [(length(xy)+1:length(xy)+length(boundary_y))' Nearest_Idx_b_y];

D = RBF_PHS_FD_all(boundary_y,xy1,Nearest_Idx_b_y,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary,cfg.rbf.derivative_order);
Dy_b_0 = D{2};

Dy_b = Dy_b_0;
Dy_b_1 = diag(Dy_b(:,Nearest_Idx_b_y(:,1)));
Dy_b(:,Nearest_Idx_b_y(:,1)) = zeros(length(boundary_y),length(boundary_y));

[Nearest_Idx_b_out] = nearest_interp(boundary_out,xy,cfg.rbf.stencil_size_boundary_cylinder);
 Nearest_Idx_b_out = [(length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out)))' Nearest_Idx_b_out];

D =  RBF_PHS_FD_all(boundary_out,xy1,Nearest_Idx_b_out,cfg.rbf.order_boundary,cfg.rbf.poly_degree_boundary, cfg.rbf.derivative_order);
Dx_b_0 = D{1};

Dx_b = Dx_b_0;
Dx_b_1 = diag(Dx_b(:,Nearest_Idx_b_out(:,1)));
Dx_b(:,Nearest_Idx_b_out(:,1)) = zeros(length(boundary_out),length(boundary_out));

D_all = RBF_PHS_FD_all(xy1,xy1,Nearest_Idx,cfg.rbf.order_main,cfg.rbf.poly_degree_main,cfg.rbf.laplacian_order);

Dx = D_all{1};
Dy = D_all{2};
L0 = D_all{3};

Nearest_Idx_nb =  xy(:,1).^2+xy(:,2).^2<(radius+r_dist)^2 ;
Idx = (1:length(xy))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);

xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),cfg.rbf.order_near_cylinder,cfg.rbf.order_near_cylinder_poly,cfg.rbf.poly_degree_near_cylinder);

Dx(Nearest_Idx_nb,:)  =  D_all_nb{1};
Dy(Nearest_Idx_nb,:)  =  D_all_nb{2};
L0(Nearest_Idx_nb,:)  =  D_all_nb{3};

% low-order PHS-RBFs and polynomials near boundaries
Nearest_Idx_nb =  xy(:,1)>x_max-x_max_dist |xy(:,1)<x_min+x_min_dist*2| xy(:,2)> y_max-y_max_dist  | xy(:,2)<y_min+y_min_dist ;
Idx = (1:length(xy))';
Nearest_Idx_nb = Idx(Nearest_Idx_nb);
xy_nb = xy1(Nearest_Idx_nb,:);

D_all_nb = RBF_PHS_FD_all(xy_nb,xy1,Nearest_Idx(Nearest_Idx_nb,:),cfg.rbf.order_near_boundary,cfg.rbf.poly_degree_main,cfg.rbf.poly_degree_near_boundary);

Dx(Nearest_Idx_nb,:)  =  D_all_nb{1};
Dy(Nearest_Idx_nb,:)  =  D_all_nb{2};
L0(Nearest_Idx_nb,:)  =  D_all_nb{3};

%% Setup time-stepping parameters and operators
nu = cfg.simulation.viscosity;  % Kinematic viscosity (1/Reynolds number)
dt = cfg.simulation.time_step;  % Time step size

% Use reduced time steps for CI environment to keep tests fast
Nt = cfg.simulation.num_time_steps;
if isCI
    Nt = cfg.simulation.num_time_steps_ci; % Reduced steps for fast CI testing
end

% Setup Laplacian operator for velocity diffusion
% Zero out Laplacian on boundary nodes (boundary conditions applied separately)
L = L0;
L(length(xy)+1:end,:) = zeros(length(boundary),length(xy1));

% Create implicit diffusion operator for Crank-Nicolson scheme
% I - (dt*nu/2)*del^2 for implicit half of viscous terms
L_I = speye(length(xy1))-dt*nu*L*cfg.schemes.crank_nicolson;

% Create separate operators for u and v velocity components
L_u = L_I;  % Copy base diffusion operator for u-velocity
L_v = L_I;  % Copy base diffusion operator for v-velocity

% Apply boundary conditions for u-velocity
L_u((length(xy)+1:length(xy)+length(boundary_y)),:) = Dy_b_0;  % Wall BC: du/dy = 0
L_u((length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out))),:) = Dx_b_0;  % Outlet BC: du/dx = 0
% Apply boundary conditions for v-velocity  
L_v((length(xy)+length(boundary_y)+1:(length(xy)+length(boundary_y)+length(boundary_out))),:) = Dx_b_0;  % Outlet BC: dv/dx = 0

% Note: Wall and inlet boundaries have Dirichlet conditions (v=0) applied directly
% Cylinder boundary has no-slip condition (u=v=0) applied directly

clear L_I L

% Precompute LU decompositions for efficient velocity solves
[LL,UU,pp,qq,rr] = lu(L_u);
L_u_inv = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));  % u-velocity solver

[LL,UU,pp,qq,rr] = lu(L_v);
L_v_inv = @(v) (qq*(UU\(LL\(pp*(rr\(v))))));  % v-velocity solver

clear L_u L_v

%% Initialize simulation variables
% Initial velocity field: uniform flow (u=1, v=0) everywhere except cylinder
V0 = zeros(length(xy1),1);  % Initial v-velocity (zero everywhere)
U0 = ones(length(xy1),1);   % Initial u-velocity (unit flow)
U0(end-length(boundary_c)+1:end) = zeros(length(boundary_c),1);  % No-slip on cylinder
W0 = [U0; V0];  % Combined velocity vector [U; V]

% Storage for velocity history (needed for multi-step time integration)
W = zeros(length(xy1)*2,Nt+1);
W(:,1) = W0;  % Store initial condition

% Initial pressure field (zero everywhere)
p0 = zeros(length(xy1_s),1);

% Define useful index lengths for boundary handling
L_B = length(boundary_c)+length(boundary_in);   % Total special boundaries
L_B_y = length(boundary_y);   % Wall boundaries
L_W = length(xy);             % Interior velocity nodes
L_B_c = length(boundary_c);   % Cylinder boundary nodes

% Define pressure boundary indices for easy access
Idx_y_s = length(xy_s)+1: length(xy_s)+length(boundary_y_s);    % Wall pressure indices
Idx_in_s = length(xy1_s)+1-length(boundary_c_s)-length(boundary_in_s): length(xy1_s)-length(boundary_c_s);  % Inlet pressure indices
Idx_out_s = length(xy_s)+1+length(boundary_y_s): length(xy_s)+length(boundary_out_s)+length(boundary_y_s);  % Outlet pressure indices

%% Main time-stepping loop: Fractional Step Method
for j = 1:Nt
    disp(['Time step j = ' num2str(j)])
    
    % Use different schemes for startup vs main integration
    if j < 3
        % First few steps: use simple first-order scheme for stability
        W(:,j+1) = W(:,j);  % Copy previous solution for startup
    else
        % Main fractional step algorithm implemented in NS_2d_cylinder_PHS
        % This performs: 1) Advection-diffusion step with Adams-Bashforth + Crank-Nicolson
        %                2) Pressure correction step to enforce incompressibility
        %                3) Velocity correction using pressure gradient
        [W(:,j+1),p0] = NS_2d_cylinder_PHS(dt,nu,W(:,j-1),W(:,j),Dy,Dx,L_inv_s,L_u_inv,L_v_inv,L0,L_B,L_B_c,L_W,L_B_y,length(boundary_s),D0_12_x,D0_12_y,D0_21_x,D0_21_y,Dy_b,Dy_b_1,D0_21_x_c,D0_21_y_c,p0,W(:,1));
    end
    % Check for numerical instability
    if isnan(W(1,j+1))
        warning('Simulation became unstable (NaN detected). Stopping at time step %d', j);
        break;
    end
end

% Extract final velocity field for potential continuation
W0 = W(:,end);

%% Visualization of final results (if plotting enabled)
if doPlot
    figure('Name','1/Re = 1e-2');
    colormap(jet)

    % Extract final time step velocity components  
    j = Nt;
    U = W(1:length(xy1),(j-1)*1+1);              % Final u-velocity
    V = W(length(xy1)+1:end,(j-1)*1+1);          % Final v-velocity
    
    % Plot u-velocity field (perturbation from uniform flow)
    subplot(2,1,1);
    % Plot u-1 to show deviation from uniform flow (u=1)
    scatter(x1,y1,cfg.visualization.scatter_size*ones(length(xy1),1),1*(U-1),'.');
    axis equal, axis tight, hold on;
    xlim([x_min x_max]);
    ylim([y_min y_max]);
    yticks(cfg.visualization.plot_tick_y)
    xticks(cfg.visualization.plot_tick_x)
    title('u-velocity perturbation (u-1)');
    ylabel('y');
    xlabel('x');
    shading interp;
    caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

    % Plot v-velocity field (should be zero in uniform flow)  
    subplot(2,1,2);
    scatter(x1,y1,cfg.visualization.scatter_size*ones(length(xy1),1),1*V,'.');
    axis equal, axis tight, hold on;
    shading interp;

    xlim([x_min x_max]);
    ylim([y_min y_max]);
    yticks(cfg.visualization.plot_tick_y)
    xticks(cfg.visualization.plot_tick_x)
    title('v-velocity');

    xlabel('x');
    ylabel('y');
    set(gca,'Ytick',[]);  % Remove y-tick labels for cleaner appearance
    caxis([-cfg.visualization.color_axis_range cfg.visualization.color_axis_range]);

    drawnow;  % Update display immediately
end